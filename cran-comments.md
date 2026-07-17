## R CMD check results

0 errors | 0 warnings | 0 notes (`devtools::check(remote = TRUE, manual = TRUE, incoming = TRUE)`,
which enables the same incoming checks CRAN's own submission pipeline runs)

* This is a new release.

* A plain `devtools::check()` run (without `incoming = TRUE`) surfaces one WARNING
  (`'qpdf' is needed for checks on size reduction of PDFs`) and one NOTE (`unable to verify current
  time`). Both are artifacts of the local checking environment (no `qpdf` binary installed; no
  outbound access to a time server), not present under the stricter incoming-checks run above, and
  not expected to reproduce on CRAN's own check machines.
