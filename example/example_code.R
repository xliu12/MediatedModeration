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
devtools::load_all("../multiMedMod")

## alternatively, install from github
# devtools::install_github(repo = "xliu12/MediatedModeration", subdir = "multiMedMod", build_vignettes = TRUE)
# library(multiMedMod)

# ?multiMedMod::MedMod



# Load data ----

load("example_data.RData")


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


# Plot ----

library(ggpattern)


plotdf_Rfoc <- out %>%
  filter(
    str_detect(Estimand, "theta")
  ) %>% 
  mutate(
    numtype = case_when(
      Estimand %in% c(glue("theta(t{c(0,1)},r1,rjo1)"), glue("theta(t{c(0,1)},r0,rjo0)")) ~ 1 ,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,r1,r1)")) ~ 2,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,r1,r0)")) ~ 3 ,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,r0,r0)")) ~ 4 ,
      Estimand %in% c(glue("theta(t{c(0,1)},r1,rjo0)")) ~ 5 
    ),
    type = factor(numtype, levels = c(1:5), labels = c(
      "As natural", #1
      "If each of M1 and M2 were independently set at their distributions among the focal subgroup (within strata of the covariates)", #2
      "If M2 were shifted to equal the reference subgroup, while fixing M1 at its distribution among the focal subgroup (within strata of the covariates)", #3
      "If M1 were also shifted to equal the reference subgroup, while fixing M2 at its distribution among the reference subgroup (within strata of the covariates)",
      "If both M1 and M2 jointly were shifted to equal the reference subgroup (within strata of the covariates)"
    )),
    subgroup = case_when(
      substr(Estimand, 10,11) == "r1" ~ "females (focal)",
      substr(Estimand, 10,11) == "r0" ~ "males (reference)"
    ),
    `Treatment assignment` = factor(ifelse(str_detect(Estimand, "t1"), "Under condition: Perceived teacher discrimination = 1", "Under condition: Perceived teacher discrimination = 0"))
  )



fig_foc <- plotdf_Rfoc %>% 
  mutate(
    `Treatment assignment` = factor(ifelse(str_detect(Estimand, "t1"), "Perceived teacher discrimination = 1", "Perceived teacher discrimination = 0")) 
  ) %>% 
  ggplot(aes(y = Estimate, x = subgroup, 
             pattern = type #, color = type 
             , fill = type
  ) ) + 
  geom_bar_pattern(
    aes(linetype = type, fill = type),
    # fill = "white", 
    # pattern_fill = "black",
    pattern_angle = 45,
    # pattern_density = 0.05,
    pattern_spacing = 0.1,
    pattern_key_scale_factor = 0.3,
    alpha = 0.4, stat = "identity", 
    color = "black",
    linewidth =0.4,
    position = position_dodge()) +
  geom_errorbar( aes(ymin = `95% CI.lower`, ymax = `95% CI.upper`
                     ,linetype = type
  ),
  position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9), 
  width = 0.4, linewidth =0.6 ) +
  scale_y_continuous("Estimated mean score of depressive symptoms (Y)") +
  coord_cartesian(ylim=c(7,12)) + 
  scale_x_discrete("Subgroup") + 
  scale_fill_discrete("Mediators (M1, M2)") +
  scale_color_discrete("Mediators (M1, M2)") +
  scale_linetype_discrete("Mediators (M1, M2)") +
  scale_pattern_discrete("Mediators (M1, M2)")+
  facet_grid(. ~ `Treatment assignment`, labeller = label_value) +
  theme_bw() +
  theme(panel.grid.minor = element_line(linewidth = 0),
        panel.grid.major.x = element_line(linewidth = 0),
        panel.grid.major.y = element_line(linewidth = 0.5, lineend = "round", color = "grey", linetype = "longdash"),
        plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 13),
        legend.text = element_text(size = 21),
        legend.title = element_text(size = 20),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.position = "bottom",
        legend.spacing.x = unit(0.2, "mm"),
        legend.key.height = unit(10, "mm"),
        legend.key.width = unit(10, "mm"),
        legend.key.size = unit(10, "mm")) 

fig_foc

