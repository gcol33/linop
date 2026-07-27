## Which base generics can S7 attach a method to? One at a time.
library(S7)
cat("S7 version:", as.character(packageVersion("S7")), "\n")
cat("R  version:", R.version.string, "\n\n")

op <- new_class("op7probe", properties = list(m = class_any))

gens <- list(
  `%*%`       = quote(`%*%`),
  crossprod   = quote(crossprod),
  tcrossprod  = quote(tcrossprod),
  t           = quote(t),
  Conj        = quote(Conj),
  solve       = quote(solve),
  dim         = quote(dim),
  nrow        = quote(nrow),
  ncol        = quote(ncol),
  as.matrix   = quote(as.matrix),
  print       = quote(print),
  summary     = quote(summary),
  `+`         = quote(`+`),
  `*`         = quote(`*`)
)

res <- data.frame(generic = character(), ok = logical(), note = character(),
                  stringsAsFactors = FALSE)

for (nm in names(gens)) {
  g <- eval(gens[[nm]])
  ok <- TRUE; note <- ""
  tryCatch({
    if (nm %in% c("t", "Conj", "dim", "nrow", "ncol", "as.matrix", "print", "summary")) {
      method(g, op) <- function(x, ...) "dispatched"
    } else {
      method(g, list(op, class_any)) <- function(x, y, ...) "dispatched"
    }
  }, error = function(e) { ok <<- FALSE; note <<- conditionMessage(e) })
  res <- rbind(res, data.frame(generic = nm, ok = ok,
                               note = substr(gsub("[\r\n]+", " ", note), 1, 90),
                               stringsAsFactors = FALSE))
}
print(res, right = FALSE)

cat("\n--- why: what does R think these are? ---\n")
for (nm in c("%*%", "crossprod", "tcrossprod", "t", "Conj", "solve", "dim")) {
  g <- get(nm, envir = baseenv())
  cat(sprintf("%-11s primitive=%-5s  ", nm, is.primitive(g)))
  bod <- if (is.primitive(g)) "<primitive>" else paste(deparse(body(g))[1], collapse = "")
  cat("body1=", substr(bod, 1, 46), "\n", sep = "")
}

cat("\n--- S7's own test for S3 genericity ---\n")
f <- tryCatch(S7:::is_S3_generic, error = function(e) NULL)
if (!is.null(f)) print(f) else cat("S7:::is_S3_generic not found\n")

cat("\n--- does base dispatch actually work on these primitives (S3, no S7)? ---\n")
zz <- structure(list(m = matrix(1:4, 2)), class = "zzprobe")
`%*%.zzprobe`      <- function(x, y) "S3 matmul OK"
crossprod.zzprobe  <- function(x, y = NULL) "S3 crossprod OK"
tcrossprod.zzprobe <- function(x, y = NULL) "S3 tcrossprod OK"
Conj.zzprobe       <- function(z) "S3 Conj OK"
for (e in c("zz %*% zz", "crossprod(zz)", "tcrossprod(zz)", "Conj(zz)")) {
  cat(sprintf("%-18s -> %s\n", e,
      tryCatch(paste(eval(parse(text = e))), error = function(err) paste("ERR:", conditionMessage(err)))))
}
