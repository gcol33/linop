.onLoad <- function(libname, pkgname) {
  register_leaf_nodes()
  register_composite_nodes()
  register_jacobi_node()
  register_section_node()
  register_structured_nodes()
  register_fourier_nodes()
  register_chebyshev_node()
  invisible(NULL)
}
