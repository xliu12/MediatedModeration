# code for the example
library(tidyverse)
library(labelled)
library(purrr)
library(origami)
library(mvtnorm)
library(parallel)
library(tidyverse)
library(glue)
library(SuperLearner)


# download the package and install
devtools::load_all("multiMedMod")

## alternatively, install from github
# devtools::install_github(repo = "xliu12/MediatedModeration", subdir = "multiMedMod", build_vignettes = TRUE)
# library(multiMedMod)

# ?multiMedMod::MedMod



# Load data ----

load("example/example_data.RData")


data_in$ID <- 1:nrow(data_in)

Rname <- "R_sex"
ttname <- "tt_tchdis"
Yname <- "Y_depress"
Cnames <- grep("^C_", colnames(data_in), value = TRUE)
M1names <- "binM1_mastery"
M2names <- "binM2_emotsupp"



# Run ----

set.seed(12345)
out <- MedMod(data_in = data_in, 
              Yname = Yname,
              M1names = M1names, M2names = M2names, ttname = ttname, Rname = Rname, Cnames = Cnames,
              learners = c("SL.glm", "SL.mean",  "SL.ranger", "SL.glmnet", "SL.earth", "SL.nnet"),
              num_folds = 4)


out


