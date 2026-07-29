# external pack scorers via in-memory assets (closed set, CRAN-safe)

# synthetic packs over a synthetic scoring panel

fake_pcbrainage_pack <- function(cpgs, seed = 42L) {
  withr::with_seed(seed, {
    coef_vec <- stats::rnorm(length(cpgs))
    impute_vec <- stats::runif(length(cpgs))
  })
  list(
    group_id = "PCBrainAge",
    cpgs = cpgs,
    coefficient_matrix = matrix(
      coef_vec,
      ncol = 1L,
      dimnames = list(NULL, "PCBrainAge")
    ),
    impute = impute_vec
  )
}

fake_pcclocks_pack <- function(cpgs, seed = 3L) {
  members <- mc_index$clock_id[mc_index$group_id == "PCClocks"]
  withr::with_seed(seed, {
    M <- matrix(
      stats::rnorm(length(cpgs) * length(members)),
      length(cpgs),
      length(members),
      dimnames = list(NULL, members)
    )
    impute_vec <- stats::runif(length(cpgs))
  })
  list(
    group_id = "PCClocks",
    cpgs = cpgs,
    coefficient_matrix = M,
    impute = impute_vec
  )
}

# systemsAge full pack layout for the family orchestrator
fake_systemsage_pack <- function(cpgs, seed = 1L) {
  order <- unname(systemsage_stack_map("SystemsAge")) # 12 labels, stack order
  organs <- setdiff(order, "Age_prediction") # 11 organ labels
  ncpg <- length(cpgs)
  pcs <- paste0("PC", seq_len(12L))
  comp_file <- function(name) {
    comp <- Filter(
      function(x) identical(x$name, name),
      clock_components("SystemsAge")
    )
    comp[[1]]$file
  }
  withr::with_seed(seed, {
    organs_mat <- matrix(
      stats::rnorm(ncpg * 11L),
      ncpg,
      11L,
      dimnames = list(NULL, organs)
    )
    systems_mat <- matrix(
      stats::rnorm(ncpg * 11L),
      ncpg,
      11L,
      dimnames = list(NULL, organs)
    )
    age_vec <- stats::rnorm(ncpg)
    impute_vec <- stats::runif(ncpg)
    rot <- matrix(stats::rnorm(144L), 12L, 12L)
    center <- stats::setNames(stats::rnorm(12L), order)
    scale <- stats::setNames(stats::runif(12L, 0.5, 1.5), order)
    model <- stats::setNames(stats::rnorm(12L), pcs)
  })
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
fake_panel <- function(n) sprintf("cg%08d", seq_len(n))

pcba_cpgs <- fake_panel(400L)
pcba_pack <- fake_pcbrainage_pack(pcba_cpgs)

pcc_cpgs <- fake_panel(300L) # shared PCClocks panel
pcc_pack <- fake_pcclocks_pack(pcc_cpgs)
pcc_members <- mc_index$clock_id[mc_index$group_id == "PCClocks"]

sa_cpgs <- fake_panel(200L)
sa_pack <- fake_systemsage_pack(sa_cpgs)
sa_members <- mc_index$clock_id[mc_index$group_id == "SystemsAge"]

# PCBrainAge

test_that("calc_clocks() scores PCBrainAge end-to-end from an in-memory pack (closed set)", {
  DNAm <- random_betas(pcba_cpgs, n = 3L)
  res <- calc_clocks(DNAm, "PCBrainAge", ext_data = pcba_pack)
  expect_setequal(colnames(res$scores), "PCBrainAge")
  expect_equal(nrow(res$scores), 3L)
  expect_false(anyNA(res$scores))

  expect_equal(res$coverage$per_clock$PCBrainAge$score_imputed_full, 0L)
})

test_that("calc_clocks() vendor-fills absent external CpGs from the pack $impute vector", {
  drop <- pcba_cpgs[1:5]
  present <- setdiff(pcba_cpgs, drop)
  DNAm <- random_betas(pcba_cpgs, n = 3L)[, present, drop = FALSE]

  res <- calc_clocks(DNAm, "PCBrainAge", ext_data = pcba_pack)
  expect_false(anyNA(res$scores))

  expect_equal(res$coverage$per_clock$PCBrainAge$score_imputed_full, 5L)
})

test_that("calc_clocks() on an external clock errors (closed set) when its pack is absent", {
  DNAm <- random_betas(pcba_cpgs, n = 2L)

  wrong <- list(
    group_id = "PCClocks",
    cpgs = "cg0001",
    coefficient_matrix = matrix(1, 1, dimnames = list(NULL, "PCADM")),
    impute = 0
  )
  expect_error(calc_clocks(DNAm, "PCBrainAge", ext_data = wrong))
})

# PCClocks

test_that("calc_clocks('PCClocks') batches all members end-to-end (closed set)", {
  DNAm <- random_betas(pcc_cpgs, n = 4L)
  pheno <- data.frame(
    ID = rownames(DNAm),
    Age = c(40, 55, 63, 71),
    Female = c(1L, 0L, 1L, 0L)
  )
  res <- calc_clocks(DNAm, "PCClocks", pheno = pheno, ext_data = pcc_pack)
  expect_setequal(colnames(res$scores), pcc_members)
  expect_equal(nrow(res$scores), 4L)
  expect_false(anyNA(res$scores))
})

test_that("requesting a subset of PCClocks returns only those columns (no expansion)", {
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

  expect_equal(
    sub$scores[, c("PCHorvath1", "PCADM")],
    full$scores[, c("PCHorvath1", "PCADM")]
  )
})

# systemsAge

test_that("calc_clocks('SystemsAge') scores the whole group (13 cols) end-to-end (closed set)", {
  DNAm <- random_betas(sa_cpgs, n = 3L)
  res <- calc_clocks(DNAm, "SystemsAge", ext_data = sa_pack)
  expect_setequal(colnames(res$scores), sa_members)
  expect_equal(nrow(res$scores), 3L)
  expect_false(anyNA(res$scores))
})

test_that("calc_clocks() vendor-fills absent SystemsAge CpGs from the pack $impute vector", {
  drop <- sa_cpgs[1:4]
  present <- setdiff(sa_cpgs, drop)
  DNAm <- random_betas(sa_cpgs, n = 3L)[, present, drop = FALSE]

  res <- calc_clocks(DNAm, "Age_prediction", ext_data = sa_pack)
  expect_false(anyNA(res$scores))
  cov <- res$coverage$per_clock$Age_prediction
  expect_equal(cov$score_imputed_full, 4L)
})

# stack labels are declared (default and columns override)

test_that("a stack that declares no columns labels each column by its operand", {
  for (id in c("GrimAgeV1", "GrimAgeV2", "DNAmPhysAge")) {
    step <- stack_step(id)
    expect_equal(stack_labels(step, id), stack_operands(step))
  }
})

test_that("SystemsAge takes its stack labels from the declared columns", {
  step <- stack_step("SystemsAge")
  expect_equal(
    stack_labels(step, "SystemsAge"),
    as.character(unlist(step[["columns"]]))
  )
  # the labels are the group's member clocks, Age_prediction included
  expect_setequal(
    stack_labels(step, "SystemsAge"),
    setdiff(sa_members, "SystemsAge")
  )
  # ...and each one pairs with the operand at its own position
  map <- stack_label_map(step, "SystemsAge")
  expect_equal(names(map), stack_operands(step))
  expect_equal(unname(map[["ap_scaled"]]), "Age_prediction")
})

test_that("a columns list that does not cover every operand is a hard stop", {
  step <- stack_step("SystemsAge")
  step[["columns"]] <- step[["columns"]][1:3]
  expect_error(stack_labels(step, "SystemsAge"))
})

# accessors over a pack (shares the pack builders above).

test_that("external accessors read the impute vector from the pack", {
  expect_equal(
    clock_impute_ref("PCADM", list(PCClocks = pcc_pack)),
    stats::setNames(pcc_pack$impute, pcc_cpgs)
  )
})

test_that("external accessors error without the group's pack", {
  expect_error(clock_impute_ref("PCADM", list()))
  expect_error(clock_impute_ref("PCADM", NULL))
})
