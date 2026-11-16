## Build the pkgdown site into docs/.
##
## pkgdown renders every Markdown file in the package root as a page. The local
## files listed in `private` are not part of the site, so they are moved aside
## for the build and restored afterwards, which also keeps them out of the
## search index.
##
## Run from the package root: Rscript pkgdown/build-site.R

private <- c("CLAUDE.md", "implementation_plan.md")
private <- private[file.exists(private)]

hold <- file.path(dirname(normalizePath(".")), ".linop-site-hold")
dir.create(hold, showWarnings = FALSE)
moved <- file.rename(private, file.path(hold, private))
if (!all(moved)) stop("could not move aside: ", paste(private[!moved], collapse = ", "))

tryCatch(
  pkgdown::build_site(".", preview = FALSE, install = TRUE),
  finally = {
    back <- file.rename(file.path(hold, private), private)
    if (!all(back)) {
      stop("restore failed; the files are in ", hold)
    }
    unlink(hold, recursive = TRUE)
  }
)

stale <- file.path("docs", c(outer(tools::file_path_sans_ext(private),
                                   c(".html", ".md"), paste0)))
unlink(stale[file.exists(stale)])
