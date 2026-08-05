# fixture builders shared by more than one test file.

# pheno frame over `ids`. a NULL column is left out entirely.
mc_pheno <- function(ids, Age = NULL, Female = NULL) {
  ph <- data.frame(ID = ids, stringsAsFactors = FALSE)
  if (!is.null(Female)) {
    ph[["Female"]] <- as.integer(Female)
  }
  if (!is.null(Age)) {
    ph[["Age"]] <- Age
  }
  ph
}

# the usual stand-in for a real Age column
mc_ages <- function(n) seq(45, 70, length.out = n)

mc_fake_cpgs <- function(n) sprintf("cg%08d", seq_len(n))

# DNAmGrip_wAge: sex-routed alias. matrix is the union of its members' panels.
grip_fixture <- function(female = c(1, 1, 0, 0), age = NULL) {
  n <- length(female)
  DNAm <- random_betas(clock_cpgs("DNAmGrip_wAge"), n = n)
  age <- age %||% mc_ages(n)
  list(
    DNAm = DNAm,
    age = age,
    female = female,
    pheno = mc_pheno(rownames(DNAm), Age = age, Female = female),
    fem = clock_coefs("DNAmGrip_wAge_Female"),
    mal = clock_coefs("DNAmGrip_wAge_Male")
  )
}

# DNAmGait_noAge_Female with `n_drop` vendor-fillable CpGs held out.
gait_holed_fixture <- function(n_drop = 5L, n = 4L) {
  id <- "DNAmGait_noAge_Female"
  coef <- clock_coefs(id)
  medians <- clock_impute_ref(id)
  drop <- intersect(names(coef), names(medians))[seq_len(n_drop)]
  DNAm <- random_betas(clock_cpgs("DNAmGait_noAge"), n = n)
  DNAm <- DNAm[, setdiff(colnames(DNAm), drop), drop = FALSE]
  list(
    id = id,
    DNAm = DNAm,
    drop = drop,
    coef = coef,
    medians = medians,
    pheno = mc_pheno(rownames(DNAm), Age = mc_ages(n), Female = rep(1L, n))
  )
}
