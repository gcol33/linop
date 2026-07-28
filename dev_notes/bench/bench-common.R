## Shared machinery for the benchmark harness: the two units it reports, the
## operators that count what a solve asks of them, and where results land.
##
## Two units, in this order.
##
##   applies   how many times a run asked the operator for A X or A^H X, and how
##             many columns went with each. This is the quantity a matrix-free
##             caller pays for, and it is the same number on every machine.
##   seconds   wall clock, the median of BENCH_REPS runs after one untimed
##             warm-up. A fact about the machine in environment.txt and about no
##             other machine, which is why it is reported second.
##
## Applies do not price everything, and no table here should be read as if they
## did. GMRES orthogonalises every step against a stored basis of up to `restart`
## vectors, and the eigensolvers against a whole round's worth; both spend
## arithmetic that no apply count sees. Reporting both columns is what keeps that
## visible.
##
## Where a second route to the same answer exists -- a dense factorisation, a
## sparse one, base R's eigen() -- it is run and reported in the same table on
## the same rows. Several of those routes come out ahead. Printing them is the
## point: the claim a matrix-free package makes is about which operators can be
## reached at all, not about winning a factorisation it also offers to call.

BENCH_DIR <- file.path("dev_notes", "bench")
BENCH_RESULTS <- file.path(BENCH_DIR, "results")

## Timing repetitions. Raising it narrows the median and changes no other number
## in the results, and it is recorded in environment.txt so a committed table
## always says how many runs produced it.
BENCH_REPS <- as.integer(Sys.getenv("LINOP_BENCH_REPS", "3"))

bench_init <- function() {
  dir.create(BENCH_RESULTS, showWarnings = FALSE, recursive = TRUE)
  invisible(BENCH_RESULTS)
}

## The machine, committed next to the numbers it produced. A timing without one
## of these is not a measurement of anything.
bench_environment <- function() {
  soft <- tryCatch(extSoftVersion(), error = function(e) character(0))
  ## An empty path is what R reports for a build whose BLAS it does not identify,
  ## and saying so is the only honest reading: it does not name which library the
  ## dense rows below went through.
  unreported <- function(x, by) {
    if (is.na(x) || !nzchar(x)) sprintf("no path reported by %s", by) else x
  }
  get1 <- function(x, nm) if (nm %in% names(x)) unname(x[[nm]]) else NA_character_
  data.frame(
    field = c("date", "R", "platform", "os", "cpu", "BLAS", "LAPACK",
              "LAPACK version", "timing repetitions"),
    value = c(format(Sys.Date()),
              R.version.string,
              R.version$platform,
              paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
              nzchar_or_na(Sys.getenv("PROCESSOR_IDENTIFIER")),
              unreported(get1(soft, "BLAS"), "extSoftVersion()"),
              unreported(tryCatch(La_library(), error = function(e) NA_character_),
                         "La_library()"),
              tryCatch(La_version(), error = function(e) NA_character_),
              as.character(BENCH_REPS)),
    stringsAsFactors = FALSE)
}

nzchar_or_na <- function(x) if (nzchar(x)) x else NA_character_

## ------------------------------------------------------ counting an operator

## A counter is an environment so the closures below can write to it from inside
## a solve without the operator holding any state a linop is not allowed to hold.
new_counter <- function() {
  e <- new.env(parent = emptyenv())
  e$applies <- 0L
  e$cols <- 0L
  e
}

counter_reset <- function(cn) {
  cn$applies <- 0L
  cn$cols <- 0L
  invisible(cn)
}

count_calls <- function(f, cn) {
  force(f); force(cn)
  function(X) {
    cn$applies <- cn$applies + 1L
    cn$cols <- cn$cols + ncol(X)
    f(X)
  }
}

## The counted operators declare their properties rather than establishing them
## by computation, because a callback leaf is what an operator with a counter in
## it has to be. So `method = "auto"` is not what these tables measure, and every
## row names its method. What auto would choose, and the one branch of it nothing
## reaches, is dev_notes/solve-dispatch-and-the-empty-branch.md.
counted_linop <- function(apply, adjoint = NULL, dim, counter,
                          properties = NULL, dtype = "double") {
  linop(count_calls(apply, counter),
        adjoint = if (is.null(adjoint)) NULL else count_calls(adjoint, counter),
        dim = dim, dtype = dtype, properties = properties)
}

## A stored matrix reached through the same counting callback, so a dense fixture
## and a matrix-free one are measured by one instrument. The apply is the BLAS
## call it would be anywhere; what the wrapper adds is the count.
counted_dense <- function(M, counter, properties = NULL) {
  force(M)
  cplx <- is.complex(M)
  counted_linop(function(X) M %*% X,
                function(X) crossprod(Conj(M), X),
                dim = dim(M), counter = counter, properties = properties,
                dtype = if (cplx) "complex" else "double")
}

## ------------------------------------------------------------------- timing

## One untimed warm-up run, whose value and whose apply counts are the ones
## returned, then BENCH_REPS timed runs. The counter is read after the warm-up so
## a count is always one run's worth however many the timing took.
bench_run <- function(expr, counter = NULL, reps = BENCH_REPS) {
  e <- substitute(expr)
  env <- parent.frame()
  if (!is.null(counter)) counter_reset(counter)
  value <- eval(e, env)
  applies <- if (is.null(counter)) NA_integer_ else counter$applies
  cols <- if (is.null(counter)) NA_integer_ else counter$cols
  secs <- vapply(seq_len(reps), function(i) {
    t0 <- Sys.time()
    eval(e, env)
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }, numeric(1))
  list(value = value, seconds = stats::median(secs),
       applies = applies, cols = cols)
}

## Memory is reported as object.size() of what a route holds and nothing else.
## gc()'s "max used" counts allocation between collections rather than what is
## live at any moment: it reads 112 Mb for a CG solve at n = 1000 whose whole
## working set is a handful of 8 kb vectors, which is a churn figure wearing a
## peak-memory label. object.size() of the operator and of the matrix are exact
## and are the two quantities the size ladder is actually about.

## ------------------------------------------------------------------ reading

vnorm <- function(x) sqrt(sum(Mod(as.vector(x))^2))

## Relative error against a truth, in the 2-norm, over whatever shape both are.
rel_error <- function(x, truth) {
  den <- vnorm(truth)
  vnorm(as.vector(x) - as.vector(truth)) / if (den > 0) den else 1
}

## A certificate line's status, without reaching into the package for it: the
## checks table is part of the object a caller gets back.
cert_line <- function(cert, check) {
  s <- cert$checks$status[cert$checks$check == check]
  if (length(s)) s else NA_character_
}

## Which of the two backward-error readings the certificate used. Decided inside
## solve_certificate() by measurement rather than by the method that ran, so it
## is read back out the same way.
cert_reading <- function(cert) {
  d <- cert$checks$detail[cert$checks$check == "backward error"]
  if (!length(d)) return(NA_character_)
  if (grepl("A^H r", d, fixed = TRUE)) "least squares" else "compatible"
}

## ------------------------------------------------------------------- output

bench_rows <- function() {
  rows <- list()
  list(
    add = function(...) {
      rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
      invisible(NULL)
    },
    collect = function() {
      if (!length(rows)) return(NULL)
      do.call(rbind, rows)
    })
}

bench_write <- function(df, name) {
  bench_init()
  path <- file.path(BENCH_RESULTS, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows)\n", path, nrow(df)))
  invisible(path)
}

fmt_cell <- function(v) {
  vapply(seq_along(v), function(i) {
    x <- v[[i]]
    if (length(x) != 1L || is.na(x)) return("-")
    if (is.numeric(x)) {
      if (is.finite(x) && x == floor(x) && abs(x) < 1e6) {
        return(format(x, trim = TRUE, scientific = FALSE))
      }
      return(formatC(x, format = "g", digits = 3))
    }
    as.character(x)
  }, character(1))
}

bench_md_table <- function(df) {
  cells <- lapply(df, fmt_cell)
  head <- paste0("| ", paste(names(df), collapse = " | "), " |")
  rule <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- vapply(seq_len(nrow(df)), function(i) {
    paste0("| ", paste(vapply(cells, `[`, character(1), i), collapse = " | "), " |")
  }, character(1))
  c(head, rule, body)
}
