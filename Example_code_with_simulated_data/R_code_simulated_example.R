
#### R Code for the Mediated-Moderation Analysis ####

# install and load the R package -------------
devtools::install_github("xliu12/MediatedModeration", subdir = "MedMod")
library(MedMod)

# install dependency packages that are not already installed
pkgs <- c("tidyverse", "glue", "origami", "mvtnorm", "SuperLearner", "ranger", "nnet")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(tidyverse)
library(glue)
library(origami)
library(mvtnorm)
library(SuperLearner)
library(ranger)
library(nnet)

# Import data. Edit the file path from which the data are to be read from
data <- read.csv("Simulation_Demo/Example_code_with_simulated_data/simulated_dataset.csv")


# Run ----------
set.seed(12345)
out <- MedMod::MedMod(
  data = data, # data containing all variables
  outcome = "Outcome", # name of outcome
  mediators = c("M1", "M2"), # name of mediators. When more than two mediators are included, the first mediator is designated as M1, and all remaining mediators are collectively treated as M2.
  treatment = "Intervention", # name of treatment variable
  subgroup = "Male", # name of moderator subgroup variable
  covariates = c("C.1", "C.2", "C.3"), # names of baseline covariates
  learners = c("SL.mean", "SL.glm", "SL.ranger", "SL.nnet"), # methods from the SuperLearner package to include in the super learner ensemble (intercept-only model, generalized linear model, random forest, neural network).
  # To see all available methods, run: SuperLearner::listWrappers()
  num_folds = 4, # number of folds for cross-fitting
  ci.level = 0.95 # default: 95% confidence interval  
  )


# checking the positivity assumption for the treatment and subgroup, and for the binary mediators
positivity_treat_subgroup <- MedMod::MedMod_overlap(out)
positivity_mediator <- MedMod::MedMod_mediator_positivity(out)

# overlap plots for checking the positivity
positivity_treat_subgroup$plot 
positivity_mediator$plot

# summaries of estimates of the conditional probabilities
positivity_treat_subgroup$summary 
positivity_mediator$summary


# Extract results for mediated moderation
out %>%
  filter(Estimand %in% c("TotMod", "MedMod", "RemainMod", # Total moderation, Mediated moderation, Remaining  moderation
                         "MedMod_M1", "MedMod_M2", "MedMod_mu")) # When more than one mediators, the output also include: Mediated moderation via M1, Mediated moderation via M2, Mediated moderation due to mediators' mutual dependence

# Plot ----
# Note: The code is for illustration only

## Relabel estimates of expected outcomes for plotting
plotdf <- out %>%
  filter(str_detect(Estimand, "theta")) %>% # extract the outcomes compared in the mediated moderation analysis
  # recode for plotting. note: `numtype` 2, 3, and 4 only exist with more than one mediator.
  mutate(
    numtype = case_when(
      Estimand %in% c(glue("theta(t{c(0,1)},r1,rjo1)"), glue("theta(t{c(0,1)},r0,rjo0)")) ~ 1 ,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,r1,r1)")) ~ 2,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,r1,r0)")) ~ 3 ,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,r0,r0)")) ~ 4 ,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,rjo0)")) ~ 5
    ),
    type = factor(numtype, levels = c(1:5), labels = c(
      "Existing condition", #1
      "Adapted condition: Mediators M1 and M2 were independent given baseline covariates", #2
      "Adapted condition: M2 were matched to the reference subgroup with similar baseline covariates, while M1 were that of the focal subgroup", #3
      "Adapted condition: M1 and M2 were independently matched to the reference subgroup with similar baseline covariates",
      "Adapted condition: Mediators were jointly matched to the reference subgroup with similar baseline covariates"
    )),
    # separate outcomes between subgroups
    subgroup = case_when(
      substr(Estimand, 10,11) == "r1" ~ "boys (focal)",
      substr(Estimand, 10,11) == "r0" ~ "girls (reference)"
    ),
    # separate outcomes between treatment conditions
    `Treatment assignment` = factor(ifelse(str_detect(Estimand, "t1"), "Intervention Condition", "Control Condition"))
  )


## plotting outcomes considered in the mediated moderation (MedMod) and remaining moderation (RemainMod)
plotdf %>%
  filter(numtype %in% c(1, 5)) %>%
  ggplot(aes(y = Estimate, x = subgroup, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymax = CI.upper, ymin = CI.lower, linetype = type), position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9), width = 0.4, linewidth =0.6 ) +
  scale_y_continuous("Estimated expectation of the outcome") +
  scale_x_discrete("Subgroup") +
  scale_fill_manual("Mediator", values = scales::hue_pal()(5)[c(1,2)]) +
  scale_linetype_discrete("Mediator") + # scales::show_col(scales::hue_pal()(5)) # R default plotting colors
  scale_color_manual("Mediator", values = scales::hue_pal()(5)[c(1,2)]) +
  scale_linetype_discrete("Mediator") +
  facet_grid(. ~ `Treatment assignment`, labeller = label_value) +
  theme_bw() +
  theme(panel.grid.minor = element_line(linewidth = 0),
        panel.grid.major.x = element_line(linewidth = 0),
        panel.grid.major.y = element_line(linewidth = 0.5, lineend = "round", color = "grey", linetype = "longdash"),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 13),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.spacing.x = unit(0.2, "mm"),
        legend.key.height = unit(10, "mm"),
        legend.key.width = unit(10, "mm"),
        legend.key.size = unit(10, "mm"))



## plotting outcomes considered in the mediated moderation via each mediator (MedMod_M1, MedMod_M2)
plotdf %>%
  filter(numtype %in% c(1:4)) %>%
  ggplot(aes(y = Estimate, x = subgroup, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymax = CI.upper, ymin = CI.lower, linetype = type), position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9), width = 0.4, linewidth =0.6) +
  scale_y_continuous("Estimated expectation of the outcome") +
  scale_x_discrete("Subgroup") +
  scale_color_manual("Mediators (M1, M2)", values = scales::hue_pal()(5)[c(1,3,4,5)]) +
  scale_fill_manual("Mediators (M1, M2)", values = scales::hue_pal()(5)[c(1,3,4,5)]) +
  scale_linetype_discrete("Mediators (M1, M2)") +
  facet_grid(. ~ `Treatment assignment`, labeller = label_value) +
  theme_bw() +
  theme(panel.grid.minor = element_line(linewidth = 0),
        panel.grid.major.x = element_line(linewidth = 0),
        panel.grid.major.y = element_line(linewidth = 0.5, lineend = "round", color = "grey", linetype = "longdash"),
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 13),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.spacing.x = unit(0.2, "mm"),
        legend.key.height = unit(10, "mm"),
        legend.key.width = unit(10, "mm"),
        legend.key.size = unit(10, "mm"))


