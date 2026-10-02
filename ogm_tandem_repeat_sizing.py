#!/usr/bin/env python3
"""Estimate tandem-repeat size (number of repeat units) from Bionano OGM alignments (.xmap) and the query cmap (_q.cmap)."""

import re
import logging
import argparse
from pathlib import Path

import pandas as pd

DEFAULT_REF_CONTIG_ID = 9
DEFAULT_LABEL_LEFT = 6238
DEFAULT_LABEL_RIGHT = 6239
DEFAULT_REF_BASELINE_BP = 3788
DEFAULT_REPEAT_UNIT_BP = 6

CMAP_CHUNKSIZE = 10**6


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s")


def read_header(path):
    """Return the column names from the '#h' line of a Bionano file (.xmap / .cmap)."""
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#h"):
                return line[2:].strip().split("\t")
            if not line.startswith("#"):
                break
    raise ValueError(f"Header line '#h' not found in {path}")


def read_bionano_table(path, **kwargs):
    """Read an .xmap / .cmap file using the column names from its header."""
    header = read_header(path)
    return pd.read_csv(path, sep="\t", comment="#", header=None,
                       names=header, **kwargs)


def extract_flanking_labels(xmap, ref_contig_id, label_left, label_right):
    """For each molecule aligned to the target contig, extract the query SiteIDs matching the two flanking reference labels."""
    pattern = re.compile(
        rf"(?<=\(){label_left},([0-9]+)\)\({label_right},([0-9]+)(?=\))"
    )

    on_contig = xmap[xmap["RefContigID"] == ref_contig_id].copy()
    matches = on_contig["Alignment"].apply(pattern.findall)

    n_multi = int((matches.apply(len) > 1).sum())
    if n_multi:
        logging.warning("%d molecule(s) with several %d/%d matches; "
                        "only the first one is kept.",
                        n_multi, label_left, label_right)

    mask = matches.apply(len) >= 1
    on_contig = on_contig[mask].copy()
    first = matches[mask].apply(lambda m: m[0])
    on_contig["Qlabel_left"] = first.apply(lambda t: int(t[0]))
    on_contig["Qlabel_right"] = first.apply(lambda t: int(t[1]))

    return on_contig.melt(
        id_vars="QryContigID",
        value_vars=["Qlabel_left", "Qlabel_right"],
        var_name="side", value_name="SiteID",
    )


def add_positions(labels_long, cmap_path):
    """Add the position (bp) of each SiteID from the query cmap, read in chunks."""
    header = read_header(cmap_path)
    keep = ["CMapId", "SiteID", "Position"]
    pieces = []
    reader = pd.read_csv(cmap_path, sep="\t", comment="#", header=None,
                         names=header, chunksize=CMAP_CHUNKSIZE)
    for chunk in reader:
        pieces.append(
            labels_long.merge(
                chunk[keep], how="inner",
                left_on=["QryContigID", "SiteID"],
                right_on=["CMapId", "SiteID"],
            )
        )
    return pd.concat(pieces, ignore_index=True)


def compute_repeats(positioned, baseline_bp, repeat_unit_bp):
    """Compute the inter-label distance per molecule and derive the number of repeat units."""
    positioned = positioned.sort_values(["QryContigID", "Position"])
    distance = positioned.groupby("QryContigID")["Position"].diff()
    positioned = positioned.assign(
        distance_bp=distance,
        repeats=(distance - baseline_bp) / repeat_unit_bp,
    )
    result = positioned.dropna(subset=["distance_bp"])
    return result[["QryContigID", "distance_bp", "repeats"]]


def parse_args():
    p = argparse.ArgumentParser(
        description=("Repeat sizing from Bionano OGM data. "
                     "Defaults target the C9orf72 G4C2 expansion; "
                     "the target options allow any other locus to be explored."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io = p.add_argument_group("files")
    io.add_argument("--xmap", required=True, type=Path, help=".xmap alignment file")
    io.add_argument("--cmap", required=True, type=Path,
                    help="query cmap (_q.cmap)")
    io.add_argument("--output", required=True, type=Path, help="output CSV")

    tgt = p.add_argument_group("target (defaults = C9orf72 / G4C2)")
    tgt.add_argument("--ref-contig", type=int, default=DEFAULT_REF_CONTIG_ID,
                     help="reference contig ID carrying the region")
    tgt.add_argument("--label-left", type=int, default=DEFAULT_LABEL_LEFT,
                     help="left flanking reference label")
    tgt.add_argument("--label-right", type=int, default=DEFAULT_LABEL_RIGHT,
                     help="right flanking reference label")
    tgt.add_argument("--baseline-bp", type=float, default=DEFAULT_REF_BASELINE_BP,
                     help="reference inter-label distance without expansion (bp)")
    tgt.add_argument("--repeat-unit", type=int, default=DEFAULT_REPEAT_UNIT_BP,
                     help="length of the repeated motif (bp)")

    return p.parse_args()


def main():
    configure_logging()
    args = parse_args()
    logging.info("Starting")
    logging.info("Target: contig %d | labels %d/%d | baseline %g bp | motif %d bp",
                 args.ref_contig, args.label_left, args.label_right,
                 args.baseline_bp, args.repeat_unit)

    xmap = read_bionano_table(args.xmap)
    labels_long = extract_flanking_labels(
        xmap, args.ref_contig, args.label_left, args.label_right)
    logging.info("%d molecule(s) spanning the target region",
                 labels_long["QryContigID"].nunique())

    positioned = add_positions(labels_long, args.cmap)
    result = compute_repeats(positioned, args.baseline_bp, args.repeat_unit)
    logging.info("%d molecule(s) with an estimated size", len(result))

    result.to_csv(args.output, index=False)
    logging.info("Written: %s", args.output)
    logging.info("Done")


if __name__ == "__main__":
    main()
