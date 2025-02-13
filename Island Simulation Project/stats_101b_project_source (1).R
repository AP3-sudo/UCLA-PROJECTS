library(knitr)
library(tidyverse)
library(lmtest)
library(knitr)
library(ggfortify)
library(ggExtra)
library(RColorBrewer)
library(patchwork)


### === Latin Square Randomization =========

set.seed(90095)
# sample(1:5, replace = FALSE) # row order
# sample(1:5, replace = FALSE) # column order

## OUTPUT:
## 2 5 3 4 1
## 5 2 1 3 4


### === Data Import & Setup =========

proj.dat <- read.csv("project_data.csv")

LS_all <- proj.dat %>%
  rename(Order = Treatment.order) %>%
  mutate_at(c("Subject", "Order", "Treatment"), factor) %>%
  select(Subject, Order, Treatment, Time)

# subset data for Latin square models
LS1 <- LS_all[proj.dat$Latin.square.assignment == "A",]
LS2 <- LS_all[proj.dat$Latin.square.assignment == "B",]
LS3 <- LS_all[proj.dat$Latin.square.assignment == "C",]
LS4 <- LS_all[proj.dat$Latin.square.assignment == "D",]
LS5 <- LS_all[proj.dat$Latin.square.assignment == "E",]


### === Auxiliary Functions =========

plot_LS <- function(LS_data, LS_number = NULL, full.data = FALSE) {
  if (full.data) {
    title <- "Full Data"
  } else {
    title <- paste("Latin Square", LS_number)
  }
  
  plt <- ggplot(LS_data, aes(x=Time, y=Treatment, groups=Treatment)) +
    geom_vline(xintercept = mean(LS_data$Time), color = "firebrick4", lty = 2, lwd = 0.6) +
    geom_boxplot(aes(fill = Treatment), alpha = 0.8) + theme_light() +
    scale_fill_brewer(palette = "Blues") +
    stat_summary(fun=mean, geom="point", shape=20, size=2, color="firebrick", fill="firebrick") +
    theme(legend.position = "none") +
    labs(title = title, x = "Change in Run Time")
  return(plt)
}

print_anova <- function(model, kable = FALSE, title = NULL) {
  anova_df <- data.frame(anova(model))
  anova_rows <- rownames(anova_df)
  anova_output <- apply(anova_df, 2, function(col) {
    rounded <- as.character(round(col, 3))
    rounded[is.na(rounded)] <- ""
    names(rounded) <- anova_rows
    return(rounded)
  })
  
  for (i in seq(nrow(anova_output))) {
    p.val <- as.numeric(anova_output[i,5])
    if (is.na(p.val)) next
    if (p.val == 0) anova_output[i,5] <- "0.000"
    if (p.val < 0.001) anova_output[i,5] <- paste(anova_output[i,5], "**", sep = "")
    if (p.val < 0.05) anova_output[i,5] <- paste(anova_output[i,5], "*", sep = "")
  }
  
  output <- data.frame(anova_output)
  colnames(output) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  if (kable) {
    return(knitr::kable(output, row.names = TRUE, caption = title))
  } else {
    return(output)
  }
}

compute_R2 <- function(model) {
  SSE <- data.frame(anova(model))['Residuals', 'Sum.Sq']
  SST <- colSums(data.frame(anova(model)))['Sum.Sq']
  
  return(1 - SSE / SST)
}


### === Model Fitting =========

# Latin square 1
aov1 <- aov(Time ~ Treatment + Order + Subject, data = LS1)
mod.tab_1 <- model.tables(aov1)$tables$Treatment
anova.ls1 <- print_anova(aov1, kable = TRUE, title = "Latin Square 1")
diagnostics.ls1 <- autoplot(aov1, smooth.colour = "#2171B5") + theme_light()
r2_ls1 <- compute_R2(aov1)

# Latin square 2
aov2 <- aov(Time ~ Treatment + Order + Subject, data = LS2)

mod.tab_2 <- model.tables(aov2)$tables$Treatment
anova.ls2 <- print_anova(aov2, kable = TRUE, title = "Latin Square 2")
diagnostics.ls2 <- autoplot(aov2, smooth.colour = "#2171B5") + theme_light()
r2_ls2 <- compute_R2(aov2)

# Latin square 3
aov3 <- aov(Time ~ Treatment + Order + Subject, data = LS3)

mod.tab_3 <- model.tables(aov3)$tables$Treatment
anova.ls3 <- print_anova(aov3, kable = TRUE, title = "Latin Square 3")
diagnostics.ls3 <- autoplot(aov3, smooth.colour = "#2171B5") + theme_light()
r2_ls3 <- compute_R2(aov3)

# Latin square 4
aov4 <- aov(Time ~ Treatment + Order + Subject, data = LS4)

mod.tab_4 <- model.tables(aov4)$tables$Treatment
anova.ls4 <- print_anova(aov4, kable = TRUE, title = "Latin Square 4")
diagnostics.ls4 <- autoplot(aov4, smooth.colour = "#2171B5") + theme_light()
r2_ls4 <- compute_R2(aov4)

# Latin square 5
aov5 <- aov(Time ~ Treatment + Order + Subject, data = LS5)

mod.tab_5 <- model.tables(aov5)$tables$Treatment
anova.ls5 <- print_anova(aov5, kable = TRUE, title = "Latin Square 5")
diagnostics.ls5 <- autoplot(aov5, smooth.colour = "#2171B5") + theme_light()
r2_ls5 <- compute_R2(aov5)

# model with combined Latin squares
aov_full <- aov(Time ~ Treatment + Order + Subject, data = LS_all)
mod.tab_all <- model.tables(aov_full)$tables$Treatment
anova.all <- print_anova(aov_full, kable = TRUE, title = "Combined Latin Squares")
diagnostics.all <- autoplot(aov_full, smooth.colour = "#2171B5") + theme_light()
r2_all <- compute_R2(aov_full)

### --- Subsection: tests for heteroskedasticity ---

# using Breusch-Pagan test
# H0: constant variance
# H1: heteroskedasticity

bp.p_vals <- cbind(
  "LS 1" = bptest(Time ~ Treatment + Order + Subject, data = LS1)$p.value,
  "LS 2" = bptest(Time ~ Treatment + Order + Subject, data = LS2)$p.value,
  "LS 3" = bptest(Time ~ Treatment + Order + Subject, data = LS3)$p.value,
  "LS 4" = bptest(Time ~ Treatment + Order + Subject, data = LS4)$p.value,
  "LS 5" = bptest(Time ~ Treatment + Order + Subject, data = LS5)$p.value,
  "Full Model" = bptest(Time ~ Treatment + Order + Subject, data = LS_all)$p.value
)
rownames(bp.p_vals) <- "p-Value"
bp.p_vals.print <- kable(bp.p_vals, row.names = TRUE, caption = "BP Test Results", digits = 4)


### === Results Section =========

# PASTE INTO RMD FILE:
# effects.print, tukey.plot, tukey.print

## box plots for each Latin square
# (plot_LS(LS1, 1) + plot_LS(LS2, 2)) / (plot_LS(LS3, 3) + plot_LS(LS4, 4))
# (plot_LS(LS5, 5) + plot_LS(LS_all, full.data = TRUE)) / (plot_spacer() + plot_spacer())


# table of treatment effects
effects.df <- bind_rows(
  mod.tab_1,
  mod.tab_2,
  mod.tab_3,
  mod.tab_4,
  mod.tab_5,
  mod.tab_all
) %>% data.frame()

rownames(effects.df) <- c(
  "LS 1", "LS 2", "LS 3", "LS 4", "LS 5", "Full Data"
)

effects.print <- kable(effects.df, row.names = TRUE, caption = "Table of Treatment Effects")

# table of R^2 values
r2.df <- bind_cols(
  "LS 1" = r2_ls1,
  "LS 2" = r2_ls2,
  "LS 3" = r2_ls3,
  "LS 4" = r2_ls4,
  "LS 5" = r2_ls5,
  "Full Data" = r2_all,
)

r2.print <- kable(r2.df, caption = "Table of R-Squared Values", digits = 3)


# post-hoc analysis with Tukey HSD

aov_all.tukey <- data.frame(TukeyHSD(aov_full)[["Treatment"]])
tukey.plot <- aov_all.tukey %>%
  rownames_to_column("Treatment") %>%
  ggplot(aes(x=diff, y=Treatment, group = Treatment)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50", alpha = 0.7) +
  geom_errorbar(aes(xmin=lwr, xmax=upr), width=0.4, linewidth=1, color="#91BFDB") +
  geom_point(color="#014990", size=4) +
  labs(title="95% Family-Wise Confidence Level",
       x="Difference in Mean Levels of Treatment") +
  theme_light() + removeGrid(y = FALSE)

aov_all.tukey.df <- aov_all.tukey %>%
  round(digits = 4) %>%
  arrange(p.adj) %>%
  rownames_to_column("Pairings") %>%
  mutate(across(diff:p.adj, as.character)) %>%
  rename("Difference" = diff,
         "Lower Bound" = lwr,
         "Upper Bound" = upr,
         "p-Value, Adjusted" = p.adj)

tukey.print <- knitr::kable(aov_all.tukey.df)






