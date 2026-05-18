# quantmsdiann reanalysis manuscript

Source repository for the quantmsdiann reanalysis paper.

## Analyses

| Dataset    | Script                                              | Description                                                              |
| ---------- | --------------------------------------------------- | ------------------------------------------------------------------------ |
| PXD003539  | `analysis/figure_original_vs_quantmsdiann.py`       | NCI-60 PCT-SWATH (Guo 2019) — original vs quantmsdiann ID counts plot.   |

## Running an analysis

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r analysis/requirements.txt
python analysis/figure_original_vs_quantmsdiann.py
```

Outputs land in `analysis/figures/`. Downloaded inputs are cached under
`data/<PXD-id>/` and are git-ignored.
