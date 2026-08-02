.onLoad <- function(libname, pkgname) {
  register_leaf_nodes()
  register_composite_nodes()
  register_jacobi_node()
  register_section_node()
  invisible(NULL)
}
