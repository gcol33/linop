## Fair re-test: correct S7 signatures for solve / Ops. Then find the root cause
## of the crossprod rejection.
library(S7)
op <- new_class("op7probe2", properties = list(m = class_any))

try_m <- function(label, expr) {
  ok <- TRUE; note <- ""
  tryCatch(eval(expr), error = function(e) { ok <<- FALSE; note <<- conditionMessage(e) })
  cat(sprintf("%-34s %-5s %s\n", label, ok, substr(gsub("[\r\n]+", " ", note), 1, 80)))
  invisible(ok)
}

cat("--- fair signatures ---\n")
try_m("solve, single dispatch",  quote({ method(solve, op) <- function(a, b, ...) "ok" }))
try_m("+ with (e1, e2)",         quote({ method(`+`, list(op, class_any)) <- function(e1, e2) "ok" }))
try_m("* with (e1, e2)",         quote({ method(`*`, list(op, class_any)) <- function(e1, e2) "ok" }))
try_m("crossprod, single disp",  quote({ method(crossprod, op) <- function(x, y = NULL) "ok" }))
try_m("crossprod via S7 wrapper",quote({ method(S7::new_generic("crossprod", "x"), op) <- function(x, ...) "ok" }))

cat("\n--- does dispatch actually FIRE for the ones that registered? ---\n")
a <- op(m = matrix(1:4, 2))
for (e in c("solve(a)", "a + 1", "a * 2", "a %*% a")) {
  cat(sprintf("%-12s -> %s\n", e,
      tryCatch(paste(eval(parse(text = e))), error = function(err) paste("ERR:", substr(conditionMessage(err),1,70)))))
}

cat("\n--- root cause: S7's allowlist of base generics ---\n")
nms <- ls(envir = asNamespace("S7"), all.names = TRUE)
cand <- grep("base|primitive|generic|known|ops|math", nms, value = TRUE, ignore.case = TRUE)
cat("candidate internals:", paste(cand, collapse = ", "), "\n\n")
for (n in cand) {
  v <- get(n, envir = asNamespace("S7"))
  if (is.character(v) && length(v) > 3) {
    cat("###", n, "( length", length(v), ")\n")
    print(v)
    cat("   contains crossprod? ", "crossprod" %in% v,
        " | contains %*%? ", "%*%" %in% v, "\n\n", sep = "")
  }
}
