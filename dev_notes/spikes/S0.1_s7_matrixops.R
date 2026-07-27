library(S7)
ns <- asNamespace("S7")

cat("=== S7 internals: what does it know about matrix ops? ===\n")
for (n in c("base_matrix_ops", "internal_generics", "group_generics", "base_ops")) {
  v <- get(n, envir = ns)
  cat("###", n, ": ", class(v), " len=", length(v), "\n", sep = "")
  if (is.function(v)) print(v) else print(utils::head(v, 30))
  cat("\n")
}

cat("=== matrixOps.S7_object ===\n")
print(get("matrixOps.S7_object", envir = ns))

cat("\n=== does crossprod reach S7 dispatch via the group generic? ===\n")
op <- new_class("mo", properties = list(m = class_any))
a <- op(m = matrix(c(1, 2, 3, 4), 2))

# route 1: explicit %*% methods (known to register)
method(`%*%`, list(op, class_any)) <- function(x, y) "S7 %*% fired"
cat("a %*% 1        ->",
    tryCatch(paste(a %*% 1), error = function(e) paste("ERR:", substr(conditionMessage(e), 1, 60))), "\n")

# route 2: can we register on the matrixOps group generic itself?
r2 <- tryCatch({ method(matrixOps, list(op, class_any)) <- function(x, y) "S7 matrixOps fired"; "registered" },
               error = function(e) paste("ERR:", substr(conditionMessage(e), 1, 70)))
cat("register matrixOps ->", r2, "\n")

# route 3: plain S3 method on the S7 class name, sidestepping S7's registrar
cn <- S7::S7_class(a)@name
cat("S7 class name:", cn, "| class(a):", paste(class(a), collapse = ", "), "\n")
assign(paste0("crossprod.", cn), function(x, y = NULL) "S3-on-S7-class fired", envir = globalenv())
registerS3method("crossprod", cn, get(paste0("crossprod.", cn), envir = globalenv()),
                 envir = globalenv())
cat("crossprod(a)   ->",
    tryCatch(paste(crossprod(a)), error = function(e) paste("ERR:", substr(conditionMessage(e), 1, 60))), "\n")

cat("\n=== same three routes, plain S3 class (control) ===\n")
z <- structure(list(m = matrix(c(1, 2, 3, 4), 2)), class = "zz3")
matrixOps.zz3 <- function(x, y) paste("S3 group generic fired:", .Generic)
registerS3method("matrixOps", "zz3", matrixOps.zz3, envir = globalenv())
for (e in c("z %*% z", "crossprod(z)", "tcrossprod(z)", "crossprod(z, z)")) {
  cat(sprintf("%-18s -> %s\n", e,
      tryCatch(paste(eval(parse(text = e))), error = function(err) paste("ERR:", substr(conditionMessage(err), 1, 60)))))
}
