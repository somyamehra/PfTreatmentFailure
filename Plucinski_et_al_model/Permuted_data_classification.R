library(dplyr)
library(tidyr)
library(readr)
library(gtools)
library(readxl) 
library(parallel)

setwd(this.path::here())

source("Plucinski_et_al_code.R")

if (!dir.exists("Permuted_datasets")) dir.create("Permuted_datasets")

N_PERMUTATIONS <- 500
N_RUNS_MCMC <- 10000
N_CORES <- detectCores()

set.seed("26092024")

Angola_data <- bind_rows(read_xlsx("../Angola_TES_data/Dimbu_Angola_TES_genotypes.xlsx", skip=3), 
                         read_xlsx("../Angola_TES_data/Dimbu_Angola_TES_genotypes.xlsx", sheet=2, skip=3)) %>%
  extract(`Sample.ID`, c("Code", "Timepoint"), "([A-Z]{2}[0-9-]{6,})(D[0-9]+)", remove = F) %>%
  mutate(Timepoint=ifelse(Timepoint=="D0", "Day 0", "Day Failure"),
         Sample.Num=1:nrow(.)) %>%
  select(-Plucinski_posterior) %>% as.data.frame

isolates <- split(Angola_data[, c("Sample.ID", "Code", "Timepoint", "Sample.Num")], 
                  Angola_data$Site) %>% lapply(function(x) {rownames(x) <- x$Sample.ID; split(x, x$Timepoint)})

generate_derangement <- function(i) {
  lapply(names(isolates), function(site) {
    x <- isolates[[site]]
    is_derangement <- TRUE
    while (is_derangement) {
      my_permutation <- sample(x[["Day 0"]]$Sample.ID, nrow(x[["Day Failure"]]))
      is_derangement <- any(x[["Day Failure"]]$Code==x[["Day 0"]][my_permutation, ]$Code)
    }
    x[["Day Failure"]]$New_ID <- paste0(site, "_Paired_", 1:nrow(x[["Day Failure"]]), " Day Failure")
    x[["Day 0"]]$New_ID <- paste0(site, "_Additional_", 1:nrow(x[["Day 0"]]), " Day 0")
    x[["Day 0"]][my_permutation, ]$New_ID <- paste0(site, "_Paired_", 1:nrow(x[["Day Failure"]]), " Day 0")
    permuted_IDs <- x %>% bind_rows %>% 
      transmute(Type=ifelse(grepl("Additional_", New_ID), "Additional", "Paired"), 
                Sample.ID=New_ID, Sample.Num=Sample.Num)
    return(permuted_IDs)}) %>% bind_rows() %>% arrange(Sample.Num) %>% 
    transmute(Type=Type, Sample.ID=Sample.ID)
}

if (!file.exists("Permuted_datasets/isolate_derangements.rds")) {
  isolate_derangements <- lapply(1:N_PERMUTATIONS, generate_derangement)
  write_rds(isolate_derangements, "Permuted_datasets/isolate_derangements.rds", compress = "gz")
} else {
  isolate_derangements <- read_rds("Permuted_datasets/isolate_derangements.rds")
}

Angola_data <- Angola_data %>% select(-Sample.Num)

record_Plucinski_classification <- function(i) {
  if (!file.exists(paste0("Permuted_datasets/Plucinski_classifier_perm_", i, ".rds"))) {
    y <- run_Plucinski_classifier(bind_cols(isolate_derangements[[i]], 
                                            Angola_data[,-c(1:3)]) %>% arrange(Sample.ID), N_RUNS_MCMC)
    write_rds(y, paste0("Permuted_datasets/Plucinski_classifier_perm_", i, ".rds"), compress = "gz")
  }
}

mclapply(1:N_PERMUTATIONS, record_Plucinski_classification,
         mc.cores = N_CORES, mc.preschedule = TRUE, mc.set.seed = TRUE)

