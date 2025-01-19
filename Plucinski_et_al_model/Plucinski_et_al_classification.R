library(dplyr)
library(tidyr)
library(readr)
library(gtools)
library(readxl)
library(parallel)

setwd(this.path::here())

source("Plucinski_et_al_code.R")

if (!dir.exists("Full_dataset")) dir.create("Full_dataset")

N_RUNS_MCMC <- 100000

set.seed("26092024")

Angola_data <- 
  bind_rows(bind_cols(Type="Paired", read_xlsx("../Angola_TES_data//Dimbu_Angola_TES_genotypes.xlsx", skip=3)), 
            bind_cols(Type="Additional", read_xlsx("../Angola_TES_data/Dimbu_Angola_TES_genotypes.xlsx", sheet=2, skip=3))) %>%
  extract(`Sample.ID`, c("Code", "Timepoint"), "([A-Z]{2}[0-9-]{6,})(D[0-9]+)", remove = F) %>%
  mutate(Timepoint=ifelse(Timepoint=="D0", "Day 0", "Day Failure"),
         Sample.ID=paste0(Code, " ", Timepoint)) %>%
  select(-Plucinski_posterior, -Code, -Timepoint) %>% as.data.frame

Plucinski_classification <- list()
Plucinski_classification[["7NMS"]] <-  run_Plucinski_classifier(Angola_data, N_RUNS_MCMC)
Plucinski_classification[["Exc_TA109"]] <- 
  run_Plucinski_classifier(Angola_data[, !grepl("TA109", colnames(Angola_data))], 
                           N_RUNS_MCMC, locirepeats=c(2,2,3,3,3,3))

write_rds(Plucinski_classification, "Full_dataset/Plucinski_classifications.rds")
