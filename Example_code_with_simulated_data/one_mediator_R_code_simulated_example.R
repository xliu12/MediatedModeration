
#### R Code for Running the Mediated-Moderation Analysis with Two Mediators ####
library(tidyverse)
library(glue)
library(origami)
library(mvtnorm)
library(SuperLearner)
library(ranger)
library(gam)
library(nnet)

# Import data
data <- read.csv("../Simulation_Demo/data_one_mediator.csv")


# run ----------
out <- MedMod(data = data,
              outcome = "Outcome",
              mediators = c("M"),
              treatment = "Intervention",
              subgroup = "Male",
              covariates = c("C.1", "C.2", "C.3"),
              learners = c("SL.glm"),
              num_folds = 1)


out$estimates
