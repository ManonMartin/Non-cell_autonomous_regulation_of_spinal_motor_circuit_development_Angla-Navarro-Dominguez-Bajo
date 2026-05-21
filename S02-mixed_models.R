# ==============================================
# Mixed models to explain the number of cells
# ==============================================

rm(list=ls())

library(tidyverse)
library(pander)
library(lme4)
library(lmerTest)
library(ggplot2)
library(purrr)

# Create directories
if(!dir.exists("results/")){
  dir.create("results/")
}

if(!dir.exists("figures/")){
  dir.create("figures/")
}


# Useful functions ================

summary_dat <- function(dat){
  dat |> summarize(
    n = n(),
    moyenne = mean(number_of_cells, na.rm = TRUE),
    minimum = min(number_of_cells, na.rm = TRUE),
    mediane = median(number_of_cells, na.rm = TRUE),
    maximum = max(number_of_cells, na.rm = TRUE),
    nNA = sum(is.na(number_of_cells))
  )
}


# mark="Foxp2"
# lev="Brach"

mixMods_mark_lev <- function(dat, mark, lev){
  
  df <- dat |> dplyr::filter(marker == mark & level == lev)
  
  df <- df |> dplyr::mutate(log_number_of_cells = case_when(
    sum(which(number_of_cells == 0)) > 0 ~  log(number_of_cells + 1), 
    sum(which(number_of_cells == 0)) == 0 ~  log(number_of_cells)
  ))
  
  # print(head(df))
  
  res_lmer <- lmerTest::lmer(log_number_of_cells ~ cond + (1|embryo), 
                             data = df, na.action = "na.omit",
                             control=lmerControl(optimizer="bobyqa"))
  
  
  cf <- confint(res_lmer, method="boot")
  
  cf_sig01 <- cf[".sig01",]
  
  lm_fit <- lm(log_number_of_cells ~ cond , 
               data = df)
  
  
  if (cf_sig01[1]>0){
    cat("======================================================== \n")
    cat(paste0(mark," - ", lev, ": ", "model: lmer(log(number_of_cells) ~ cond + (1|embryo)) \n"))
    mod <- res_lmer
    mod_name <- "lmer"
    ## check model hypotheses
    # plot(mod)
    # plot(predict(mod),residuals(mod))
    # qqnorm(residuals(mod))
    # library(performance)
    # performance::check_model(res_lmer)
    cat("======================================================== \n")
    
  }else{
    cat("======================================================== \n")
    cat(paste0(mark," - ", lev, ": ", "model: lm(log(number_of_cells) ~ cond) \n"))
    mod <- lm_fit
    mod_name <- "lm"
    #   # check model hypotheses
    # plot(mod)
    # plot(predict(mod),residuals(mod))
    # qqnorm(residuals(mod))
    cat("======================================================== \n")
  }
  
  print(summary(mod))
  
  
  t <- summary(mod)
  
  pval <- t$coefficients[2,"Pr(>|t|)"]
  
  test <- c(model = mod_name, pvalue = pval)
  return(test)
}

mixMods_mark <- function(dat, mark){
  
  df <- dat |> dplyr::filter(marker == mark)
  
  df <- df |> dplyr::mutate(log_number_of_cells = case_when(
    sum(which(number_of_cells == 0)) > 0 ~  log(number_of_cells + 1), 
    sum(which(number_of_cells == 0)) == 0 ~  log(number_of_cells)
  ))
  
  # print(head(df))
  res_lmer <- lmerTest::lmer(log_number_of_cells ~ cond + (1|embryo), 
                             data = df, 
                             control=lmerControl(optimizer="bobyqa"))
  
  
  cf <- confint(res_lmer, method="boot")
  
  cf_sig01 <- cf[".sig01",]
  
  lm_fit <- lm(log_number_of_cells ~ cond , 
               data = df)
  
  if (cf_sig01[1]>0){
    cat("======================================================== \n")
    cat(paste0(mark," - ", lev, ": ", "model: lmer(log(number_of_cells) ~ cond + (1|embryo)) \n"))
    mod <- res_lmer
    mod_name <- "lmer"
    ## check model hypotheses
    # plot(mod)
    # plot(predict(mod),residuals(mod))
    # qqnorm(residuals(mod))
    # library(performance)
    # performance::check_model(res_lmer)
    cat("======================================================== \n")
    
  }else{
    cat("======================================================== \n")
    cat(paste0(mark, ": ", "model: lm(log(number_of_cells) ~ cond) \n"))
    mod <- lm_fit
    mod_name <- "lm"
    #   # check model hypotheses
    # plot(mod)
    # plot(predict(mod),residuals(mod))
    # qqnorm(residuals(mod))
    cat("======================================================== \n")
  }
  
  print(summary(mod))
  
  
  t <- summary(mod)
  
  pval <- t$coefficients[2,"Pr(>|t|)"]
  
  test <- c(model = mod_name, pvalue = pval)
  return(test)
}


#  List of figures ================

# List of figures
list.files("data", recursive = FALSE)

#  Fig1B ================
fig_name <- "Fig1B"

## Data importation
dat <- readxl::read_xlsx("data/Fig.1B/Fig1v.xlsx", sheet = 1)

dat <- dat |> dplyr::select(-Mean)

writexl::write_xlsx(dat, 
                    path = paste0("results/S02", fig_name,"_combined.xlsx"))

## Data exploration

table(dat$cond, dat$level)
table(dat$level, dat$cut)
table(dat$embryo, dat$cut)

dat_sum <- dat |> 
  group_by(marker, level) |> 
  summary_dat()

dat_sum

writexl::write_xlsx(dat_sum, 
                    path = paste0("results/S02", fig_name,"_summary.xlsx"))

dat |> 
  ggplot(aes(x = number_of_cells)) + 
  geom_density() + 
  theme_bw()

dat |> 
  ggplot(aes(x = log(number_of_cells))) + 
  geom_density() + 
  theme_bw()

dat |> ggplot(aes(x = log(number_of_cells), fill = cond)) + 
  geom_density(alpha = 0.6) + 
  facet_wrap(~marker+level, scales = "free_y", ncol = 3)+ 
  theme_bw()

dat |> mutate(log_number_of_cells = log(number_of_cells)) |> 
  select(marker, level, cond, number_of_cells)

dat |>
  ggplot(aes(x = level, y = number_of_cells, fill = cond)) + 
  geom_boxplot(position=position_dodge(1)) + 
  facet_wrap(~marker, scales = "free_y", ncol = 2) + 
  theme_bw() + 
  scale_fill_manual(values = c("#ECA700","#0067AC")) +
  xlab("Condition") +
  ylab("Number of cells")


dat |>
  ggplot(aes(x = level, y = number_of_cells,color = cond)) + 
  geom_boxplot(position=position_dodge(1), outlier.shape = NA) + 
  geom_jitter(position=position_dodge(1)) +
  facet_wrap(~marker, scales = "free_y", ncol = 2) + 
  theme_bw() + 
  scale_color_manual(values = c("#ECA700","#0067AC")) +
  xlab("Condition") +
  ylab("Number of cells")

dat |>
  ggplot(aes(x = cond, y = number_of_cells, 
             color = as.factor(embryo), shape = cut)) + 
  geom_point(position=position_dodge(1)) +
  # facet_wrap(vars(marker,level), scales = "free_y") +
  facet_grid(marker~level, scales = "free_y") +
  geom_vline(xintercept = 1.5, linetype = 2, color = "grey") + 
  theme_bw() 

# Mixed models

# - mixed model if random effect is significant
# - linear regression if random effect is not significant

ddf <- dat |> select(marker, level) |> unique()

res_map2 <- map2_df(.x = ddf$marker, 
                    .y = ddf$level, 
                    .f = ~ mixMods_mark_lev(mark = .x,lev = .y, dat = dat))

res <- cbind(ddf, res_map2)

res <- res |> 
  mutate(pvalue = as.numeric(pvalue), 
         adj_pvalue_fdr = p.adjust(pvalue, method = "fdr")) |> 
  arrange(pvalue) 

res

writexl::write_xlsx(res, path = paste0("results/S02",fig_name,".xlsx"))

# Fig1A ================

fig_name <- "Fig1A"

# Data importation

dat <- readxl::read_xlsx("data/Fig1A_Foxp2-MafA OTP-Prdm8.xlsx", sheet = 1)

str(dat)

dat <- dat |> mutate(cut = as.factor(cut))

writexl::write_xlsx(dat, 
                    path = paste0("results/S02", fig_name,"_combined.xlsx"))

# Data exploration

table(dat$cond, dat$cut)
table(dat$embryo, dat$cut)
table(dat$cut)

dat$cut <- factor(dat$cut)

dat_sum <- dat |> 
  group_by(marker) |> 
  summary_dat()

dat_sum

writexl::write_xlsx(dat_sum, 
                    path = paste0("results/S02", fig_name,"_summary.xlsx"))

dat |> 
  ggplot(aes(x = number_of_cells)) + 
  geom_density() + 
  theme_bw()

dat |> 
  ggplot(aes(x = log(number_of_cells))) + 
  geom_density() + 
  theme_bw()

dat |> 
  ggplot(aes(x = log(number_of_cells), fill = cond)) + 
  geom_density(alpha = 0.6) + 
  facet_wrap(~marker, scales = "free_y", ncol = 3)+ 
  theme_bw()

dat |>
  ggplot(aes(x = cond, y = number_of_cells, fill = cond)) + 
  geom_boxplot(position=position_dodge(1)) + 
  facet_wrap(~marker, scales = "free_y", ncol = 2) + 
  theme_bw() + 
  scale_fill_manual(values = c("#ECA700","#0067AC")) +
  xlab("Condition") +
  ylab("Number of cells")


dat |>
  ggplot(aes(x = cond, y = number_of_cells,color = cond)) + 
  geom_boxplot(position=position_dodge(1), outlier.shape = NA) + 
  geom_jitter(position=position_dodge(1)) +
  facet_wrap(~marker, scales = "free_y", ncol = 2) + 
  theme_bw() + 
  scale_color_manual(values = c("#ECA700","#0067AC")) +
  xlab("Condition") +
  ylab("Number of cells")

dat |>
  ggplot(aes(x = cond, y = number_of_cells, 
             color = as.factor(embryo), shape = cut)) + 
  geom_point(position=position_dodge(1)) +
  facet_wrap(vars(marker), scales = "free_y") +
  geom_vline(xintercept = 1.5, linetype = 2, color = "grey") + 
  theme_bw() 

# Mixed models

# - mixed model if random effect is significant
# - linear regression if random effect is not significant

ddf <- dat |> select(marker) |> unique()


res_map2 <- map_df(.x = ddf$marker, 
                   .f = ~ mixMods_mark(dat = dat, mark = .x))



res <- cbind(ddf, res_map2)

res <- res |> 
  mutate(pvalue = as.numeric(pvalue), 
         adj_pvalue_fdr = p.adjust(pvalue, method = "fdr")) |> 
  arrange(pvalue) 

res

writexl::write_xlsx(res, path = paste0("results/S02",fig_name,".xlsx"))

