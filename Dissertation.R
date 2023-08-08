## Data Wrangling: Stopsack PRAD
## Lewis Quayle, PhD
## 2023-07-27


# setup -------------------------------------------------------------------

# load libraries
library("tidyverse")
library("here")
library("janitor")
library("survival")
library("survminer")

# directory hierarchy
project_dir <- here()
in_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "proc_data")


# read data ---------------------------------------------------------------

sm_424 <- read_tsv(file = file.path(in_dir, "prad_mcspc_mskcc_2020_clinical_data.tsv"))
lg_2069 <- read_delim(file = file.path(in_dir, "prad_msk_stopsack_2021_clinical_data.tsv"))


# clean and join data -----------------------------------------------------

# clean variable names
sm_424_cln <- 
  sm_424 %>% 
  janitor::clean_names()

lg_2069_cln <- 
  lg_2069 %>% 
  janitor::clean_names()

# join data sets
joined_df <- 
  sm_424_cln %>% 
  left_join(lg_2069_cln, by = "patient_id") %>% 
  mutate(race_category = coalesce(race_category.x, race_category.y)) # fills in NA by merging two race_category cols

category_counts <- table(joined_df_cln$race_category)
print(category_counts)



# there are 24 non-concordant patient IDs so this is as good as it gets - could filter later to remove missingness
joined_df %>%
  select(starts_with("race_")) %>% 
  summarise(across(everything(), ~ sum(is.na(.))))

# clean joined data frame: select vars of interest/remove duplicate vars, re-code chr to fct/lgcl etc., filter rows
joined_df_cln <- 
  joined_df %>% 
  select(everything()) %>% # replace everything() with only vars you want to keep
  mutate(
    across(where(is_character), as_factor),
    survival_status = ifelse(str_detect(survival_status, fixed("dead", ignore_case = TRUE)), yes = 0, no = 1)
    
    # could include cleaning of variables to homogenise the format here before converting all to factor
    # may need to re-code some vars that should not be fct here
    
  ) %>% 
  filter(!is.na(race_category))


# fit survival model ------------------------------------------------------

fit <- 
  survfit(
    Surv(overall_survival_since_sample_collection_months, survival_status) ~ race_category,
    data = joined_df_cln
  )

print(fit)

# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

joined_df_cln_tab <- data.frame(time = fit$time,
                                n.risk = fit$n.risk,
                                n.event = fit$n.event,
                                n.censor = fit$n.censor,
                                surv = fit$surv,
                                upper = fit$upper,
                                lower = fit$lower
)
head(joined_df_cln_tab)


# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw())


# Get unique categories
categories <- unique(joined_df_cln$race_category)
print(categories)


ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 20,     # break X axis in time intervals by 20.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("White", "Black or African American", "Other", "Black", "Asian"),    # change legend labels.

)

summary(fit)$table

#he survival curves can be shorten using the argument xlim
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           xlim = c(0, 80))

#to plot cumulative events
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           fun = "event")

#To plot cumulative hazard
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           fun = "cumhaz")

#Kaplan-Meier life table: summary of survival curves.......................
summary(fit)

#function surv_summary() [in survminer package] to get a summary of survival curves.
res.sum <- surv_summary(fit)
head(res.sum)

attr(res.sum, "table")

#Log-Rank test comparing survival curves: survdiff()
surv_diff <- survdiff(Surv(overall_survival_since_sample_collection_months, survival_status) ~ race_category,
                      data = joined_df_cln
)
surv_diff


