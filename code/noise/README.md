# TCR dynamics (private working repository)

This is a private, working repository for clonotype noise modeling
and longitudinal dynamics analysis.

The repository is designed to be:
- implementable at every stage
- runnable on lab machines and HPC clusters
- not public / not presentation-ready until paper submission

## Structure

- `code/noise/`    : noise model and noise inference pipeline
- `code/dynamics/` : clonotype trajectory construction and dynamics inference
- `methods/`       : manuscript-ready Methods sections (DOCX)
- `examples/`      : runnable scripts (lab / cluster)
- `sandbox/`       : scratch and exploratory work

## Minimal usage

```bash
python code/noise/noise_pipeline_full.py \
  --data-dir path/to/replicates \
  --results-dir results/noise \
  --null subject \
  --alpha 0.01
