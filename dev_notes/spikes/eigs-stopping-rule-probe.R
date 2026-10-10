## What the eigenpair certificate reads after a run, on every engine, for the
## requests whose stopping rule and certificate could disagree (gcol33/linop#2).
##
## The certificate holds a pair to two rows: the backward error
## ||A x - mu x|| / (||A|| ||x||) and the residual ||A x - mu x|| / (|mu| ||x||),
## each with the arithmetic floor. A stopping rule that tests only the first locks
## pairs whose |mu| is well below ||A|| and the second row then fails; a stopping
## rule that tests on (A - sigma I)^-1 bounds the residual on A by
## tol ||A - sigma I||, not by tol |mu|. Each row below is one run: its cost and
## what its certificate says.
##
##   Rscript dev_notes/spikes/eigs-stopping-rule-probe.R <label>
##
## writes dev_notes/spikes/results/eigs-stopping-rule-<label>.csv.

label <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(label)) label <- "run"

suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-linop.R")

cert_status <- linop:::cert_status

engines <- c("lanczos", if (requireNamespace("linop.primme", quietly = TRUE)) "primme")

lap60 <- laplacian_1d(60)
lap80 <- laplacian_1d(80)

eig_cases <- list(
  list(id = "lap80 k3 SA", A = lap80, args = list(k = 3, which = "smallest_algebraic")),
  list(id = "lap60 k3 sigma=1", A = lap60, args = list(k = 3, sigma = 1)),
  list(id = "lap60 k3 sigma=2", A = lap60, args = list(k = 3, sigma = 2)),
  list(id = "lap60 k3 SA", A = lap60, args = list(k = 3, which = "smallest_algebraic")),
  list(id = "lap60 k3 LA", A = lap60, args = list(k = 3, which = "largest_algebraic")),
  list(id = "lap60 k4 SA ncv24", A = lap60,
       args = list(k = 4, which = "smallest_algebraic", ncv = 24L, maxit = 300L)),
  list(id = "lap60 k4 LA ncv24", A = lap60,
       args = list(k = 4, which = "largest_algebraic", ncv = 24L, maxit = 300L)),
  list(id = "shifted lap40 k2 LM", A = shifted_laplacian_1d(40, 2),
       args = list(k = 2, which = "largest")),
  list(id = "spd40 gap k3 LA",
       A = as_spd_linop(spd_prescribed(40, c(100, 99.5, 40, 10, seq(1, 0.1, length.out = 36)),
                                       seed = 7L)),
       args = list(k = 3, which = "largest_algebraic")),
  list(id = "spd40 k3 LA ncv12",
       A = as_spd_linop(spd_prescribed(40, c(30, 12, 11.8, 5, seq(2, 0.4, length.out = 36)))),
       args = list(k = 3, which = "largest_algebraic", ncv = 12L)),
  list(id = "lap400 k6 LA", A = laplacian_1d(400),
       args = list(k = 6, which = "largest_algebraic"))
)

svd_cases <- list(
  list(id = "diff81 k3 smallest", A = diff_1d(81), args = list(k = 3, which = "smallest")),
  list(id = "diff81 k3 largest", A = diff_1d(81), args = list(k = 3, which = "largest"))
)

row_of <- function(id, verb, method, fit) {
  cert <- fit$certificate
  rel <- cert$values$residual
  data.frame(case = id, verb = verb, method = method,
             iterations = fit$iterations, nconv = fit$nconv,
             overall = cert$overall,
             residual = cert_status(cert, "residual"),
             backward = cert_status(cert, "backward error"),
             worst_residual = if (all(is.na(rel))) NA_real_ else max(rel, na.rm = TRUE),
             worst_backward = max(cert$values$backward_error))
}

rows <- list()
for (cs in eig_cases) for (m in engines) {
  fit <- do.call(eigs, c(list(cs$A), cs$args, list(method = m)))
  rows[[length(rows) + 1L]] <- row_of(cs$id, "eigs", m, fit)
}
svd_engines <- c("golub-kahan", setdiff(engines, "lanczos"))
for (cs in svd_cases) for (m in svd_engines) {
  fit <- do.call(svds, c(list(cs$A), cs$args, list(method = m)))
  rows[[length(rows) + 1L]] <- row_of(cs$id, "svds", m, fit)
}

out <- do.call(rbind, rows)
print(out, row.names = FALSE, digits = 3)
dir.create("dev_notes/spikes/results", showWarnings = FALSE)
utils::write.csv(out, sprintf("dev_notes/spikes/results/eigs-stopping-rule-%s.csv", label),
                 row.names = FALSE)
