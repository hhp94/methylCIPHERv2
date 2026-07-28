# dev-only helpers (load_all)

# run cohort parity tests
test_parity <- function(filter = "fixtures-parity", ...) {
  withr::with_envvar(
    c(MC_PARITY = "1"),
    devtools::test(filter = filter, ...)
  )
}
