STR / SSTR: SOURCE PACKAGES AND PORTABLE RUN CODE
================================================
This software-only supplement contains source R packages, execution scripts,
frozen seeds and study definitions. It contains no saved numerical results,
manuscript PDFs, report tables, installed libraries or machine-specific paths.
The scenario/task/signal CSVs define inputs, not numerical result tables.

Research engines (original names and authorship retained):
  reference22: softridge 0.4.0, the reference R implementation.
  independent192: softridge 0.3.0.1 and scaleCorrection 0.1.0,
                  with RcppArmadillo numerical kernels.
The two softridge versions are deliberately installed into separate libraries.
Original package API names are retained; manuscript names are STR and SSTR.
SSTR denotes the corrected c0 scale; SSTR-LS denotes the earlier LS scale.
The unrelated Restricted Bridge Estimation (rbridge) package is not this engine.

See RUN_GUIDE.txt for installation, checks, bounded examples, full commands,
checkpoint/resume and result verification. No installation, data download or
numerical execution occurs merely by opening this archive.

Run python3 verify_sources.py after extraction. Source archives and scientific
functions are unchanged from the complete submission bundle. Only the run
adapters make the optional real-data file an explicit, checksum-pinned download.
No public repository URL or software DOI is asserted for these research engines.
