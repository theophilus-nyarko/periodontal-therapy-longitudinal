# =============================================================================
# Evaluating the Effectiveness of Surgical and Non-Surgical Periodontal Therapy
# Over Time: A Longitudinal Analysis

#
# Note: in the dataset, the plaque index column is named `plague`.
# =============================================================================

library(dplyr)
library(ggplot2)
library(nlme)
library(lattice)
library(reshape2)

periodontal <- read.csv("data/periodontal.csv")
dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Data checks
# -----------------------------------------------------------------------------
str(periodontal)
dim(periodontal)
head(periodontal)
class(periodontal$surgical)

# Missing data
sum(is.na(periodontal))

# Number of observations at each visit
table(periodontal$visit)

# -----------------------------------------------------------------------------
# 2. Descriptive statistics (baseline visit only, to avoid repetition)
# -----------------------------------------------------------------------------
baseline <- subset(periodontal, visit == 1)

table(baseline$surgical)                     # treatment groups
table(baseline$smoking, baseline$surgical)   # smoking by treatment
table(baseline$alcohol, baseline$surgical)   # alcohol by treatment

# -----------------------------------------------------------------------------
# 3. Exploratory plots
# -----------------------------------------------------------------------------

# Figure 1: mean PPD over time by treatment group
means <- aggregate(ppd ~ visit + surgical, data = periodontal, FUN = mean)

png("figures/01_mean_ppd_by_treatment.png", width = 800, height = 450)
plot(1:8, means$ppd[means$surgical == 0], type = "b", col = "blue", lwd = 2,
     xlab = "Visit (months)", ylab = "Average PPD (mm)",
     ylim = c(3, 5), main = "")
lines(1:8, means$ppd[means$surgical == 1], type = "b", col = "red", lwd = 2)
legend("topright", legend = c("Non-surgical", "Surgical"),
       col = c("blue", "red"), lty = 1, lwd = 2)
dev.off()

# Figure 2: individual patient trajectories (spaghetti plot)
png("figures/02_individual_trajectories.png", width = 800, height = 450)
print(xyplot(ppd ~ visit | factor(surgical), group = ID, data = periodontal,
             type = "b", xlab = "Visit", ylab = "PPD (mm)",
             strip = strip.custom(factor.levels = c("Non-surgical", "Surgical"))))
dev.off()

# PPD by visit and treatment
boxplot(ppd ~ visit * surgical, data = periodontal,
        xlab = "Visit", ylab = "PPD (mm)",
        names = rep(1:8, 2), col = rep(c("blue", "red"), each = 8))

# Relationships between PPD and continuous predictors
pairs(~ ppd + age + bmi + teeth + plague + gingival + bleeding,
      data = periodontal, main = "Relationships Between PPD and Predictors")

# Distribution of baseline PPD by treatment group
hist(baseline$ppd[baseline$surgical == 0], col = "blue", border = "white",
     xlab = "PPD (mm)", main = "Non-surgical", xlim = c(4, 5.5))
hist(baseline$ppd[baseline$surgical == 1], col = "red", border = "white",
     xlab = "PPD (mm)", main = "Surgical", xlim = c(4, 5.5))

# -----------------------------------------------------------------------------
# 4. Correlation among repeated measures
# -----------------------------------------------------------------------------
wide_ppd <- dcast(periodontal, ID ~ visit, value.var = "ppd")
cor(wide_ppd[, 2:9], use = "complete.obs")

# Figure 3: ACF of residuals from an independence model.
# The pattern does not decay smoothly, so AR(1) is not appropriate and an
# unstructured covariance (corSymm) is used in the models below.
fit_ind <- gls(ppd ~ visit * surgical, data = periodontal, method = "ML")

png("figures/03_acf_residuals.png", width = 800, height = 450)
acf(residuals(fit_ind), main = "ACF of Residuals")
dev.off()

# =============================================================================
# QUESTION 1: Does the rate of change in PPD differ between treatments?
# =============================================================================

# Step 1: linear vs. quadratic time trend (GLS, unstructured covariance)
fit_linear <- gls(ppd ~ visit * surgical,
                  data = periodontal,
                  correlation = corSymm(form = ~ visit | ID),
                  method = "ML")
summary(fit_linear)

fit_quad <- gls(ppd ~ poly(visit, degree = 2) * surgical,
                data = periodontal,
                correlation = corSymm(form = ~ visit | ID),
                method = "ML")
summary(fit_quad)

anova(fit_linear, fit_quad)

# Step 2: add demographic covariates one at a time
fit_m1 <- gls(ppd ~ poly(visit, degree = 2) * surgical,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_m2 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_m3 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + race,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_m4 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + race + bmi,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_m5 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + race + bmi + hba,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_m6 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + race + bmi + alcohol,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_m7 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + race + bmi +
                smoking + alcohol,
              data = periodontal,
              correlation = corSymm(form = ~ visit | ID), method = "ML")

AIC(fit_m1, fit_m2, fit_m3, fit_m4, fit_m5, fit_m6, fit_m7)

# Compare GLS with an LME model (random intercept and random slope), using ML
fit_gls_ml <- gls(ppd ~ poly(visit, degree = 2) * surgical + age,
                  data = periodontal,
                  correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_lme_ml <- lme(ppd ~ visit * surgical + age,
                  random = ~ visit | ID,
                  data = periodontal, method = "ML")

AIC(fit_gls_ml, fit_lme_ml)

# Step 3: final model (Model 1A), estimated with REML
fit_gls <- gls(ppd ~ poly(visit, degree = 2) * surgical + age,
               data = periodontal,
               correlation = corSymm(form = ~ visit | ID), method = "REML")
summary(fit_gls)

# Model 1A diagnostics
png("figures/04_model1A_residual_histogram.png", width = 800, height = 450)
hist(residuals(fit_gls), main = "Histogram of Residuals", xlab = "Residuals")
dev.off()

png("figures/05_model1A_qq_plot.png", width = 800, height = 450)
print(qqnorm(fit_gls, abline = c(0, 1)))
dev.off()

png("figures/06_model1A_residuals_vs_fitted.png", width = 800, height = 450)
plot(fitted(fit_gls), residuals(fit_gls),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red", lwd = 2)
dev.off()

# =============================================================================
# QUESTION 2: Does treatment effectiveness differ for current smokers with
# moderate-to-heavy alcohol use? (Model 2A)
# =============================================================================
fit_gls_q2 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age +
                    surgical:smoking:alcohol,
                  data = periodontal,
                  correlation = corSymm(form = ~ visit | ID), method = "REML")
summary(fit_gls_q2)

# =============================================================================
# QUESTION 3: Do periodontal conditions improve the model beyond demographics?
# =============================================================================

# Add periodontal conditions one at a time to Model 1A (ML for comparison)
fit_gls_1 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_gls_2 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_gls_3 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague + teeth,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_gls_4 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague + teeth +
                   gingival,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_gls_5 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague + teeth +
                   gingival + bleeding,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_gls_6 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague + teeth +
                   gingival + bleeding + crp,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

fit_gls_7 <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague + teeth +
                   gingival + bleeding + crp + rf,
                 data = periodontal,
                 correlation = corSymm(form = ~ visit | ID), method = "ML")

AIC(fit_gls_1, fit_gls_2, fit_gls_3, fit_gls_4, fit_gls_5, fit_gls_6, fit_gls_7)

# Likelihood ratio test: Model 1A (demographics) nested in Model 3A (+ plaque)
gls_reducedML <- fit_gls_1
gls_fullML    <- fit_gls_2
anova(gls_reducedML, gls_fullML)

# Model 3A, estimated with REML
gls_fullREML <- gls(ppd ~ poly(visit, degree = 2) * surgical + age + plague,
                    data = periodontal,
                    correlation = corSymm(form = ~ visit | ID), method = "REML")
summary(gls_fullREML)

# Model 3A diagnostics
png("figures/07_model3A_residual_histogram.png", width = 800, height = 450)
hist(residuals(gls_fullREML), main = "Histogram of Residuals", xlab = "Residuals")
dev.off()

png("figures/08_model3A_qq_plot.png", width = 800, height = 450)
print(qqnorm(gls_fullREML, abline = c(0, 1)))
dev.off()

png("figures/09_model3A_residuals_vs_fitted.png", width = 800, height = 450)
plot(fitted(gls_fullREML), residuals(gls_fullREML),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red", lwd = 2)
dev.off()
