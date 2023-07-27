install.packages("readr")
library(readr)

survdata <- read_tsv("prad_mcspc_mskcc_2020_clinical_data.tsv")
survdata1 <- read_tsv('prad_msk_stopsack_2021_clinical_data.tsv')

head(survdata)
summary(survdata)



install.packages("dplyr")
library(dplyr)
survdataMerge <- inner_join(survdata, survdata1, by = "Patient ID")


sum(is.na(survdataMerge$`Race Category.y`))
colSums(is.na(survdataMerge))

sum(is.na(survdataMerge$`Race Category.y`))
sum(is.na(survdataMerge$`Survival Status`))

# Assuming your dataset is named "survdataMerge" and the variable you want to analyze is "Race Category.y"
unique_categories <- unique(survdataMerge$`Race Category.y`)

num_categories <- length(unique_categories)
print(num_categories)

#Using dplyr package's select function:
survdataMerge1 <- survdataMerge %>% select(`Patient ID`,`Age at Diagnosis`,`Age at Diagnosis`, `Androgen Deprivation Therapy (ADT)`, `Timing of metastases`,`Time from sample to castration resistance [months]`,`Overall Survival since Sample Collection (Months)`,`Fraction Genome Altered.y`, `Race Category.y`,`Prostate-specific antigen.y`,Smoking,`Sample Type`,`TMB (nonsynonymous).y`,`Survival Status`,`Patient's Vital Status` )
sum(is.na(survdataMerge1))
colSums(is.na(survdataMerge1))

survdataMerge1

survdataMerge1 <- survdataMerge1 %>%
  mutate(`Race Category.y` = as.numeric(factor(`Race Category.y`, levels = unique(`Race Category.y`))) - 1)

# Encode "survival status" into 0 and 1
survdataMerge1 <- survdataMerge1 %>%
  mutate(`Survival Status` = recode(`Survival Status`, "Dead" = 1, .default = 0))


install.packages(c("survival", "survminer"))
library("survival")
library("survminer")

install.packages("ggplot2")
library('ggplot2')


data_type <- class(survdataMerge1$`Race Category.y`)
print(data_type)
status_data_type <- class(survdataMerge1$`Survival Status`)
print(status_data_type)



# Assuming dataset is named "survdataMerge1" and the "Survival Status" column is named "Survival_Status"
# Convert character values to logical (TRUE/FALSE)
survdataMerge1$`Survival Status` <- survdataMerge1$`Survival Status` == "Alive"  # Change "Event" to the appropriate event indicator value

summary(survdataMerge1)

survdataMerge1$`Race Category.y` <- factor(survdataMerge1$`Race Category.y`, levels = unique(survdataMerge1$`Race Category.y`))

# We want to compute the survival probability by Race Category
surv_obj <- Surv(time =survdataMerge1$`Overall Survival since Sample Collection (Months)`, event = survdataMerge1$`Survival Status`)



survdataMerge1



survdataMerge2 <- survdataMerge1 %>% select(`Patient ID`,`Overall Survival since Sample Collection (Months)`,`Race Category.y`,Smoking,`Survival Status`,`Patient's Vital Status` )

# Fit the survival model
surv_obj <- Surv(time =survdataMerge2$`Overall Survival since Sample Collection (Months)`, event = survdataMerge2$`Survival Status`)

surv_fit <- survfit(surv_obj ~ `Race Category.y`, data = survdataMerge2)


fit <- survfit(Surv(`Overall Survival since Sample Collection (Months)`, `Survival Status`) ~ `Race Category.y`, data = survdataMerge1)


print(fit)


# Assuming your dataset is named "survdataMerge1"
print(class(survdataMerge1$`Overall Survival since Sample Collection (Months)`))
print(class(survdataMerge1$`Survival Status`))
print(class(survdataMerge1$`Race Category.y`))
# Assuming your dataset is named "survdataMerge1"
print(names(survdataMerge1))


