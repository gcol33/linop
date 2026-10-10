## Every suite except test-backends.R is about the reference engines and the
## certificates built on them, so method = "auto" resolves to the reference
## whether or not linop.primme is installed. test-backends.R names its engines,
## and tests auto's choice with this option cleared.
options(linop.backend_priority = c(primme = -Inf))
