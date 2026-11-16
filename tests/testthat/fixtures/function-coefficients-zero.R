structure(list(subject = "chebyshev solve", checks = structure(list(
    check = c("finite residual", "truncation residual", "backward error", 
    "forward error"), status = c("pass", "pass", "pass", "not_checked"
    ), source = c("computation", "theorem", "computation", "computation"
    ), guarantee = c("identity", "deterministic_bound", "identity", 
    "identity"), confidence = c(1, 1, 1, NA), detail = c("||boundary and kept rows of L c - F|| / ||F|| = 0.000e+00 on 32 coefficients", 
    "the rows the truncation dropped, computed exactly: 0.000e+00 relative; widths 32", 
    "the coefficients solve the infinite system with the right-hand side moved by 0.000e+00 relative", 
    "a distance to the true solution needs a bound on the inverse of the operator, which is not computed"
    )), row.names = c(NA, -4L), class = "data.frame"), overall = "qualified", 
    without_deterministic_bound = "forward error", probes = NULL, 
    values = list(finite = 0, truncation = 0, backward_error = 0, 
        n = 32L, widths = 32L), evidence = NULL), class = "linop_certificate")
