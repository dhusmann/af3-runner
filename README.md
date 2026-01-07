# af3-runner

High-throughput AlphaFold3 pipeline scripts for Sherlock/Slurm, plus utilities to organize and sync outputs (including Google Drive via `rclone`).

## Where to look

- Pipeline: `README-AF3PipelineForSherlock.md`
- Output organization + sync: `README-OutputSyncAndOrganization.md`
- Parallel sync deployment notes: `PARALLEL_SYNC_DEPLOYMENT.md`
- Full script index: `SCRIPT_SUMMARY.md`
- Job creation: `createAF3query.sh` (run `./createAF3query.sh --help`)

## Getting started (typical)

1. Clone/copy this repo into your **run directory** on scratch (e.g. `/scratch/groups/.../alphafold3/`).
2. Copy `examples/folding_jobs.example.csv` to `folding_jobs.csv` in the run directory and edit job names.
3. Follow `README-AF3PipelineForSherlock.md` to launch/monitor the pipeline.
