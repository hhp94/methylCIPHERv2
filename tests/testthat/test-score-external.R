# external pack scorers via in-memory assets (closed set, CRAN-safe). smoke only.

# synthetic packs over a synthetic scoring panel

fake_pcbrainage_pack <- function(cpgs) {
  list(
    group_id = "PCBrainAge",
    cpgs = cpgs,
    coefficient_matrix = matrix(
      stats::rnorm(length(cpgs)),
      ncol = 1L,
      dimnames = list(NULL, "PCBrainAge")
    ),
    impute = stats::runif(length(cpgs))
  )
}

fake_pcclocks_pack <- function(cpgs, members) {
  list(
    group_id = "PCClocks",
    cpgs = cpgs,
    coefficient_matrix = matrix(
      stats::rnorm(length(cpgs) * length(members)),
      length(cpgs),
      length(members),
      dimnames = list(NULL, members)
    ),
    impute = stats::runif(length(cpgs))
  )
}

# systemsAge full pack layout. dimensions come from the declared stack.
fake_systemsage_pack <- function(cpgs) {
  order <- unname(systemsage_stack_map("SystemsAge")) # labels, stack order
  organs <- setdiff(order, "Age_prediction")
  ncpg <- length(cpgs)
  n_sys <- length(order)
  n_org <- length(organs)
  pcs <- paste0("PC", seq_along(order))
  comp_file <- function(name) {
    comp <- Filter(
      function(x) identical(x$name, name),
      clock_components("SystemsAge")
    )
    comp[[1]]$file
  }
  organ_mat <- function() {
    matrix(
      stats::rnorm(ncpg * n_org),
      ncpg,
      n_org,
      dimnames = list(NULL, organs)
    )
  }
  organs_mat <- organ_mat()
  systems_mat <- organ_mat()
  age_vec <- stats::rnorm(ncpg)
  impute_vec <- stats::runif(ncpg)
  rot <- matrix(stats::rnorm(n_sys^2L), n_sys, n_sys)
  center <- stats::setNames(stats::rnorm(n_sys), order)
  scale <- stats::setNames(stats::runif(n_sys, 0.5, 1.5), order)
  model <- stats::setNames(stats::rnorm(n_sys), pcs)
  rot_df <- cbind(
    data.frame(system = order, stringsAsFactors = FALSE),
    stats::setNames(as.data.frame(rot), pcs)
  )
  list(
    group_id = "SystemsAge",
    cpgs = cpgs,
    organs = organs_mat,
    systems = systems_mat,
    age = age_vec,
    impute = impute_vec,
    tensors = stats::setNames(
      list(center, scale, rot_df, model),
      c(
        comp_file("systems_pca_center"),
        comp_file("systems_pca_scale"),
        comp_file("systems_pca_rotation"),
        comp_file("systems_model")
      )
    )
  )
}

# synthetic pack panels (closed set, sizes differ to catch cross-wiring)
pcba_cpgs <- mc_fake_cpgs(400L)
pcba_pack <- fake_pcbrainage_pack(pcba_cpgs)

pcc_cpgs <- mc_fake_cpgs(300L) # shared PCClocks panel
pcc_members <- mc_index$clock_id[mc_index$group_id == "PCClocks"]
pcc_pack <- fake_pcclocks_pack(pcc_cpgs, pcc_members)

sa_cpgs <- mc_fake_cpgs(200L)
sa_pack <- fake_systemsage_pack(sa_cpgs)
sa_members <- mc_index$clock_id[mc_index$group_id == "SystemsAge"]

# pcBrainAge

test_that("calc_clocks() scores PCBrainAge end-to-end from an in-memory pack (closed set)", {
  skip_on_cran()
  DNAm <- random_betas(pcba_cpgs, n = 3L)
  res <- calc_clocks(DNAm, "PCBrainAge", ext_data = pcba_pack)
  expect_setequal(colnames(res$scores), "PCBrainAge")
  expect_equal(nrow(res$scores), 3L)
  expect_false(anyNA(res$scores))

  expect_equal(res$coverage$per_clock[[1]]$PCBrainAge$score_imputed_full, 0L)
})

test_that("calc_clocks() vendor-fills absent external CpGs from the pack $impute vector", {
  skip_on_cran()
  drop <- pcba_cpgs[1:5]
  present <- setdiff(pcba_cpgs, drop)
  DNAm <- random_betas(pcba_cpgs, n = 3L)[, present, drop = FALSE]

  res <- calc_clocks(DNAm, "PCBrainAge", ext_data = pcba_pack)
  expect_false(anyNA(res$scores))

  expect_equal(res$coverage$per_clock[[1]]$PCBrainAge$score_imputed_full, 5L)
})

test_that("calc_clocks() refuses a supplied pack belonging to another group", {
  skip_on_cran()
  DNAm <- random_betas(pcba_cpgs, n = 2L)

  wrong <- list(
    group_id = "PCClocks",
    cpgs = "cg0001",
    coefficient_matrix = matrix(1, 1, dimnames = list(NULL, "PCADM")),
    impute = 0
  )
  expect_error(calc_clocks(DNAm, "PCBrainAge", ext_data = wrong))
})

# pcClocks

test_that("requesting a subset of PCClocks returns only those columns (no expansion)", {
  skip_on_cran()
  DNAm <- random_betas(pcc_cpgs, n = 3L)
  pheno <- data.frame(
    ID = rownames(DNAm),
    Age = c(50, 60, 70),
    Female = c(0L, 1L, 1L)
  )

  sub <- calc_clocks(
    DNAm,
    c("PCHorvath1", "PCADM"),
    pheno = pheno,
    ext_data = pcc_pack
  )
  full <- calc_clocks(DNAm, "PCClocks", pheno = pheno, ext_data = pcc_pack)

  expect_setequal(colnames(sub$scores), c("PCHorvath1", "PCADM"))
  expect_setequal(colnames(full$scores), pcc_members)
  expect_false(anyNA(full$scores))

  # the batched member score does not depend on who else was requested
  expect_equal(
    sub$scores[, c("PCHorvath1", "PCADM")],
    full$scores[, c("PCHorvath1", "PCADM")]
  )
})

# systemsAge

test_that("calc_clocks('SystemsAge') scores the whole group end-to-end (closed set)", {
  skip_on_cran()
  DNAm <- random_betas(sa_cpgs, n = 3L)
  res <- calc_clocks(DNAm, "SystemsAge", ext_data = sa_pack)
  expect_setequal(colnames(res$scores), sa_members)
  expect_equal(nrow(res$scores), 3L)
  expect_false(anyNA(res$scores))
})
