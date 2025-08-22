rm(list = ls())

library(SKM)
library(tidyverse)
library(foreach)
library(doParallel)
library(BGLR)

# GLOBAL PARAMS ----------------------------------------------------------------
dataset_file <- "EYT_1"
cv_folds_num <- 10
cores_num <- 10
cv_testing_prop <- 0.5
iterations_number <- 10000
burn_in <- 5000
base_results_dir <- "results"

# TEST PARAMS ------------------------------------------------------------------
dataset_file <- "Test"
cv_folds_num <- 2
cores_num <- 2
cv_testing_prop <- 0.5
iterations_number <- 10000
burn_in <- 5000
base_results_dir <- "trash"

decompose_matrices <- function(K, A) {
  KA <- K %*% A

  UT <- KA
  UT[!upper.tri(KA)] <- 0
  CC <- UT + t(UT)
  diag(CC) <- diag(KA)
  CC <- CC / mean(diag(CC))

  LT <- KA
  LT[!lower.tri(KA)] <- 0
  PP <- LT + t(LT)
  diag(PP) <- diag(KA)
  PP <- PP / mean(diag(PP))

  return(list(CC = CC, PP = PP))
}

# PREPARE DATA -----------------------------------------------------------------
load(sprintf("data/%s.RData", dataset_file), verbose = TRUE)

Pheno <- Pheno %>%
  mutate(GID = Line) %>%
  arrange(Env, GID) %>%
  droplevels()

# Select the unique lines both in pheno and geno and sort them
final_geno_lines <- sort(intersect(Pheno$GID, rownames(Geno)))
Geno <- Geno[final_geno_lines, final_geno_lines]

results_dir <- file.path(base_results_dir, dataset_file)

all_predictions <- list()

cluster <- makeCluster(cores_num, outfile = "")
registerDoParallel(cluster)

folds <- cv_random_line(
  Pheno$Line,
  folds_number = cv_folds_num,
  testing_proportion = cv_testing_prop
)

ZL <- model.matrix(~ 0 + Line, data = Pheno)
ZE <- model.matrix(~ 0 + Env, data = Pheno)
Geno <- data.matrix(Geno)
K_E <- ZE %*% t(ZE)
K_G <- ZL %*% Geno %*% t(ZL)
K_GE <- K_E * K_G

K_C <- decompose_matrices(K = K_E, A = K_G)
K_CC <- K_C$CC
K_PP <- K_C$PP

for (t in seq_along(data_info$traits)) {
  trait <- data_info$traits[t]
  SKM::echo(
    "*** Trait '%s': %i / %i ***",
    trait,
    t,
    length(data_info$traits)
  )

  y <- Pheno[[trait]]

  fold_predictions <- foreach(fold = folds, .combine = append) %dopar% {
    SKM::echo("\t*** Fold: %i / %i ***", fold$num, length(folds))

    `%>%` <- dplyr::`%>%`

    y_na <- replace(y, fold$testing, NA)
    observed <- y[fold$testing]

    predictions <- list()

    # ETA0 -------------------------------------------------------------------
    ETA0 <- list(
      Env = list(model = "RKHS", K = K_E),
      Lines = list(model = "RKHS", K = K_G)
    )

    model0 <- BGLR::BGLR(
      y = y_na,
      ETA = ETA0,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- append(predictions, list(data.frame(
      Model = "A",
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = observed,
      Predicted = model0$yHat[fold$testing]
    )))

    # ETA1 -------------------------------------------------------------------
    ETA1 <- list(
      Env = list(model = "RKHS", K = K_E),
      Lines = list(model = "RKHS", K = K_G),
      GE = list(model = "RKHS", K = K_GE)
    )

    model1 <- BGLR::BGLR(
      y = y_na,
      ETA = ETA1,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- append(predictions, list(data.frame(
      Model = "B",
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = observed,
      Predicted = model1$yHat[fold$testing]
    )))

    # ETA2 -------------------------------------------------------------------
    ETA2 <- list(
      Env = list(model = "RKHS", K = K_E),
      Lines = list(model = "RKHS", K = K_G),
      K_CC = list(model = "RKHS", K = K_CC)
    )

    model2 <- BGLR::BGLR(
      y = y_na,
      ETA = ETA2,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- append(predictions, list(data.frame(
      Model = "C",
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = observed,
      Predicted = model2$yHat[fold$testing]
    )))

    # ETA3 -------------------------------------------------------------------
    ETA3 <- list(
      Env = list(model = "RKHS", K = K_E),
      Lines = list(model = "RKHS", K = K_G),
      K_PP = list(model = "RKHS", K = K_PP)
    )

    model3 <- BGLR::BGLR(
      y = y_na,
      ETA = ETA3,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- append(predictions, list(data.frame(
      Model = "D",
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = observed,
      Predicted = model3$yHat[fold$testing]
    )))

    # ETA4 -------------------------------------------------------------------
    ETA4 <- list(
      Env = list(model = "RKHS", K = K_E),
      Lines = list(model = "RKHS", K = K_G),
      K_CC = list(model = "RKHS", K = K_CC),
      K_PP = list(model = "RKHS", K = K_PP)
    )

    model4 <- BGLR::BGLR(
      y = y_na,
      ETA = ETA4,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- append(predictions, list(data.frame(
      Model = "E",
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = observed,
      Predicted = model4$yHat[fold$testing]
    )))

    # ETA5 -------------------------------------------------------------------
    ETA5 <- list(
      Env = list(model = "RKHS", K = K_E),
      Lines = list(model = "RKHS", K = K_G),
      GE = list(model = "RKHS", K = K_GE),
      K_CC = list(model = "RKHS", K = K_CC),
      K_PP = list(model = "RKHS", K = K_PP)
    )

    model5 <- BGLR::BGLR(
      y = y_na,
      ETA = ETA5,
      response_type = "gaussian",
      nIter = iterations_number,
      burnIn = burn_in,
      verbose = FALSE,
      saveAt = tempdir()
    )
    predictions <- append(predictions, list(data.frame(
      Model = "F",
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = observed,
      Predicted = model5$yHat[fold$testing]
    )))

    # Merge fold predictions -------------------------------------------------
    predictions <- dplyr::bind_rows(predictions) %>%
      dplyr::mutate(
        Dataset = data_info$name,
        Trait = trait,
        Fold = fold$num
      )

    list(predictions)
  }

  all_predictions <- append(all_predictions, fold_predictions)
}

Predictions <- bind_rows(all_predictions) %>%
  round_df(4)

Summary <- Predictions %>%
  group_by(Dataset, Trait, Env, Model, Fold) %>%
  summarise(
    Cor = cor(Observed, Predicted),
    MSE = mse(Observed, Predicted),
    NRMSE = nrmse(Observed, Predicted),
    PM_10 = best_lines_match(tibble(Line, Observed, Predicted), 10),
    PM_20 = best_lines_match(tibble(Line, Observed, Predicted), 20),
    PM_30 = best_lines_match(tibble(Line, Observed, Predicted), 30)
  ) %>%
  summarise_if(
    is.numeric,
    list(
      MEAN = ~ mean(., na.rm = TRUE),
      SE = ~ sd(., na.rm = TRUE) / sqrt(n())
    ),
    .groups = "keep"
  ) %>%
  round_df(4) %>%
  select(-contains("Fold")) %>%
  as_tibble() %>%
  rename_with(~gsub("_MEAN", "", .))

mkdir(results_dir)
write_csv(Predictions, file.path(results_dir, "predictions.csv"))
write_csv(Summary, file.path(results_dir, "summary.csv"))
