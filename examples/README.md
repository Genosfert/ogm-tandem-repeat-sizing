# Example data

`example.xmap` and `example_q.cmap` are **small synthetic files** written to test
the installation. They do not come from real samples.

They contain four molecules:

| QryContigID | Content | Expected repeats |
|-------------|---------|------------------|
| 101 | spans the default target, distance = baseline | 0 |
| 102 | spans the default target, +45 units | 45 |
| 103 | spans the default target, +100 units | 100 |
| 104 | aligned elsewhere (contig 12), must be ignored | – |

Run:

```bash
python ogm_tandem_repeat_sizing.py \
    --xmap examples/example.xmap \
    --cmap examples/example_q.cmap \
    --output results.csv
```

The result must be identical to `expected_output.csv`.
