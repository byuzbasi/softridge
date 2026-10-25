# Reproduce the reference and independent studies

The `reproduction/` folder in the GitHub repository contains software
and run definitions. It has no saved performance reports. Download the
repository and run commands from its reproduction directory:

``` bash
cd reproduction
python3 verify_sources.py
Rscript --vanilla RUNME.R --help
Rscript --vanilla RUNME.R install --study reference22
Rscript --vanilla RUNME.R install --study independent192
Rscript --vanilla RUNME.R check --study reference22
Rscript --vanilla RUNME.R check --study independent192
```

The runners install archived research packages into separate project
libraries. They do not download third-party dependencies. A local source
fingerprint records which engine is used. `RUN_GUIDE.txt` lists the
declared dependencies, recorded versions and explicit
execution/resume/verification commands.

| Study | Scenarios | Replications per scenario | Engine |
|:---|---:|---:|:---|
| Reference | 22 | 100 | Original R softridge 0.4.0 |
| Independent | 192 | 500 | RcppArmadillo softridge 0.3.0.1 + scaleCorrection 0.1.0 |

The studies retain their own seeds, designs and tuning protocols. The
12-by-40 reference grid differs from the independent 20-by-41 fold-local
grid. The different version numbers identify historical implementations;
they do not imply that the pure-R reference supersedes the compiled
research engine.

Only explicit `--all --confirm-full` commands start complete studies.
The independent runner has ten validated batch checkpoints, a single
writer, signature-checked resume and checksum verification. Existing
outputs are protected. End-to-end portable full-run cost has not been
measured; inspect a bounded example before planning an independent full
reproduction.

## Data dictionary

| Input | Meaning |
|:---|:---|
| `config.rds` | Frozen numerical study definition and seed list |
| `scenarios.csv` | Sample size, dimension, sparsity, AR correlation and SNR |
| `tasks.csv` | Scenario-replication identities and recorded order |
| `signals.csv` | Active-variable locations and coefficients |
| `protocol_corrected.json` | Corrected-scale protocol and metric definitions |

The optional riboflavin example uses public hdi data. It is downloaded
only through the explicit `fetch_riboflavin.R` command and checked
against a pinned SHA-256 digest. Original data are retained without
modification.
