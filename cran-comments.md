## R CMD check results

0 errors | 1 warning | 1 note

* This is a new release.

* The WARNING (`'qpdf' is needed for checks on size reduction of PDFs`) and
  the NOTE (`unable to verify current time`) are artifacts of the local
  checking environment (no `qpdf` binary installed; no outbound access to a
  time server) and are not expected to reproduce on CRAN's own check
  machines.
