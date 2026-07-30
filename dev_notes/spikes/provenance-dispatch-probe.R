## Can a provider reach its own method through the four provenance generics as core
## stores the envelope? Section 5.11 has core carrying an opaque envelope and never
## inspecting it, and R/provenance.R exports provenance_lift(), provenance_refine(),
## provenance_original_residual() and provenance_summary() as S3 generics dispatching
## on the envelope. The envelope is built by set_provenance() as a bare
## list(provider =, payload =), so this checks what class the dispatch actually sees.
##
## Run from the repository root:
##   "C:/Program Files/R/R-4.6.0/bin/x64/Rscript.exe" dev_notes/spikes/provenance-dispatch-probe.R

devtools::load_all(".", quiet = TRUE)

out_dir <- file.path("dev_notes", "spikes", "results")
if (!dir.exists(out_dir)) stop("run from the repository root")
txt <- file(file.path(out_dir, "provenance-dispatch.txt"), open = "wt")
say <- function(...) { l <- paste0(...); cat(l, "\n", sep = ""); cat(l, "\n", sep = "", file = txt) }

A <- set_provenance(linop(diag(3)), "linop.hilbert", list(scheme = "finite_section", n = 128))
p <- provenance(A)

say("class(envelope) as stored: ", paste(class(p), collapse = ", "))
say("attributes(envelope):      ", paste(names(attributes(p)), collapse = ", "))
say("")

## A provider defines a method the obvious way and tries to reach it.
registerS3method("provenance_lift", "hilbert_provenance",
                 function(p, x, ...) "provider method reached")

r_stored <- tryCatch(provenance_lift(p, 1:3),
                     error = function(e) paste("ERROR:", conditionMessage(e)))
r_classed <- tryCatch(provenance_lift(structure(p, class = "hilbert_provenance"), 1:3),
                      error = function(e) paste("ERROR:", conditionMessage(e)))

say("provenance_lift(envelope as stored):  ", r_stored)
say("provenance_lift(envelope hand-classed): ", r_classed)
say("")

## Does a class on the payload survive set_provenance()?
A2 <- set_provenance(linop(diag(3)), "linop.hilbert",
                     structure(list(scheme = "fs"), class = "hilbert_payload"))
say("class(payload) after set_provenance:  ",
    paste(class(provenance(A2)$payload), collapse = ", "))
say("class(envelope) after set_provenance: ",
    paste(class(provenance(A2)), collapse = ", "))
say("")

## What dispatching on the payload instead would do, without changing the signature.
## The generic stands in for the one in R/provenance.R; the method is defined where
## this generic can see it, which inside the package would be registerS3method().
prov_lift_v2 <- function(p, x, ...) UseMethod("prov_lift_v2", p$payload)
prov_lift_v2.default <- function(p, x, ...) stop("no method")
prov_lift_v2.hilbert_payload <- function(p, x, ...) "provider method reached via payload class"
say("UseMethod on p$payload:               ",
    tryCatch(prov_lift_v2(provenance(A2), 1:3),
             error = function(e) paste("ERROR:", conditionMessage(e))))

close(txt)
