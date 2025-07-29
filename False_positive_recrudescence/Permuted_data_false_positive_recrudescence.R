library(readr)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(reshape2)
library(tidyr)
library(stringr)
library(ggplot2)
library(parallel)

library(showtext)
font_add(family = "Arial", regular = "Arial.ttf") ## here is the path to the font to add.
showtext.auto()

library(extrafont)
font_import() # takes a few minutes
loadfonts(device="postscript")

library(PfRecur)
library(PDQutils)
library(copula)
library(poisbinom)
library(VGAM)

setwd(this.path::here())

# ============================ PARSE PERMUTED DATA ===============================

Angola_2021_data <- read_rds("../Angola_TES_data/Dimbu_Angola_TES_data.rds")

MARKERS <- c(M313="M313", M383="M383", TA1="TA1", POLYA="POLYA",
             PFPK2="PFPK2", M2490="M2490", TA109="TA109") 

EPSILON_DEFAULT <- 0.05
OMEGA_VALS <- seq(0.75, 1, 0.05)
OMEGA_DEFAULT <- 0.9
BETA <- 0.25

N_CORES <- detectCores()

isolate_derangements <- read_rds("../Plucinski_et_al_model/Permuted_datasets/isolate_derangements.rds")

permuted_isolates <- isolate_derangements %>% 
  lapply(function(x) x %>% subset(Type=="Paired") %>% 
           extract(Sample.ID, c("New_Code", "Timepoint"), "([A-Za-z0-9_]+) (Day [A-Za-z0]+)") %>% 
           tibble::rownames_to_column(var="Sample_ID") %>% 
           dcast(New_Code~Timepoint, value.var="Sample_ID")) %>% 
  bind_rows(.id="iter") %>% mutate(New_Code=gsub("Sul", "Lunda Sul", New_Code))

match_counting <- apply(permuted_isolates, 1, function(x) {
  apply(Angola_2021_data[c(x["Day 0"], x["Day Failure"]), MARKERS], 2, function(y) {
    alleles_1 <- str_split(y[1], "/")[[1]]
    alleles_1 <- alleles_1[alleles_1!=""]
    alleles_2 <- str_split(y[2], "/")[[1]]
    alleles_2 <- alleles_2[alleles_2!=""]
    return(c(n_match=length(intersect(alleles_1, alleles_2))>=1,
             n_comparable=(length(alleles_1)*length(alleles_2))>0))}) %>% rowSums()})

permuted_isolates <- bind_cols(permuted_isolates, t(match_counting))

Plucinski_permuted <- lapply(1:length(isolate_derangements), function(i) 
  read_rds(paste0("../Plucinski_et_al_model/Permuted_datasets/Plucinski_classifier_perm_", i, ".rds")))

permuted_isolates <- lapply(Plucinski_permuted, function(x) 
  x[["posterior_recrudesence"]] %>% transmute(New_Code=isolate, Plucinski=prob)) %>% 
  bind_rows(.id="iter") %>% merge(permuted_isolates)

# ============================ PARSE DATA ===============================

set.seed("20042024")

day_0_isolates <- unique(permuted_isolates$`Day 0`)
day_0_isolates <- split(day_0_isolates, substr(day_0_isolates, 1, 1))

Angola_2021_marker_set <- lapply(MARKERS, function(m) {
  sort(unique(setdiff(unlist(strsplit(Angola_2021_data[,m], "/")), "NA")))})

Angola_2021_genotype_matrix <- lapply(MARKERS, function(m) {
  keep_samples <- subset(Angola_2021_data[, "Sample_ID"], 
                         Angola_2021_data[, m]!="")
  GT <- matrix(0, nrow=length(keep_samples), ncol=length(Angola_2021_marker_set[[m]]),
               dimnames = list(keep_samples, Angola_2021_marker_set[[m]]))
  for (indiv in keep_samples) {
    GT[indiv, unlist(setdiff(strsplit(Angola_2021_data[indiv, m], "/"), "NA"))] <- 1
  }
  return(GT)
})

Angola_2021_locus_cardinality <- lapply(Angola_2021_genotype_matrix, rowSums)


# ============================ PARSE DATA ===============================

error_matrix <- PfRecur::geometric_error_matrix(Angola_2021_marker_set, EPSILON_DEFAULT)
  
classify_pair <- function(i) {
  x <- PfRecur::evaluate_posterior(permuted_isolates[i, "Day Failure"], 
                                   permuted_isolates[i, "Day 0"],
                                   setdiff(day_0_isolates[[substr(permuted_isolates[i, "Day 0"], 1, 1)]], 
                                           permuted_isolates[i, "Day 0"]),
                                   Angola_2021_genotype_matrix,
                                   error_matrix, OMEGA_VALS, BETA, MARKERS)
  return(bind_cols(x[["metrics"]], permuted_isolates[i,]))
}

if (!file.exists("Dimbu_Angola_permuted_classifications.rds")) {
  PfRecur_permuted_data <-
    mclapply(1:nrow(permuted_isolates), classify_pair,
             mc.cores = N_CORES, mc.preschedule = TRUE, mc.set.seed = TRUE) %>% 
    bind_rows() %>% mutate(epsilon=EPSILON_DEFAULT) %>% 
    extract(New_Code, "Site", "([A-Za-z ]+)_.")
  write_rds(PfRecur_permuted_data, "Dimbu_Angola_permuted_classifications.rds", compress = "gz")
} else {
  PfRecur_permuted_data <- read_rds("Dimbu_Angola_permuted_classifications.rds")
}

match_counting_fpr <- permuted_isolates %>% subset(n_comparable==7) %>%
  extract(New_Code, "Site", "([A-Za-z ]+)_.") %>% 
  group_by(Site, iter, n_match) %>% summarise(count=n()) %>% split(f=.$Site~.$iter) %>% 
  lapply(function(x) {x <- x %>% arrange(-n_match); x$prop <- x$count/sum(x$count); 
  x$ecdf=cumsum(x$prop); return(x)}) %>% bind_rows() %>% 
  subset(n_match<=5 & n_match>3) %>%
  bind_rows(data.frame(Site="Lunda Sul", iter=unique(permuted_isolates$iter),
                       n_match=5, count=0, prop=0, ecdf=0)) %>%
  transmute(Site=Site, iter=iter, FPR=ecdf, model=paste0("≥", n_match, "/7"))

PfRecur_fpr <- PfRecur_permuted_data %>% 
  group_by(Site, omega, epsilon, iter) %>% 
  summarise(FPR=mean(M1))

PfRecur_fpr_default <- PfRecur_fpr %>% as.data.frame %>%
  subset(omega==OMEGA_DEFAULT & epsilon==EPSILON_DEFAULT) %>%
  select(Site, iter, FPR) %>% mutate(model="PfRecur")

Plucinski_fpr <- PfRecur_permuted_data %>% 
  select(Site, iter, `Day 0`, `Day Failure`, Plucinski) %>% unique %>% 
  group_by(Site, iter) %>% 
  summarise(FPR=mean(Plucinski)) %>% mutate(model="CDC") %>% as.data.frame

fpr_summary_plot <- bind_rows(PfRecur_fpr_default, match_counting_fpr, Plucinski_fpr) %>% 
  group_by(Site, model) %>% 
  summarise(median=median(FPR), 
            LCI=quantile(FPR, 0.025), UCI=quantile(FPR, 0.975),
            LCI2=quantile(FPR, 0.25), UCI2=quantile(FPR, 0.75)) %>% 
  mutate(model=factor(model, levels=c("PfRecur", "CDC", "≥4/7", "≥5/7"))) %>% 
  ggplot(aes(x=model, y=median)) + 
  geom_point(cex=2) +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.55, col="#737373") + 
  geom_errorbar(aes(ymin=LCI2, ymax=UCI2), width=0.3, lwd=0.65) + 
  facet_grid(cols=vars(Site)) +
  scale_y_continuous(breaks=seq(0, 0.3, 0.02)) +
  xlab("Classification method") +
  ylab("False positive\nrecrudescence rate") +
  theme_bw() + theme(panel.grid.major.x = element_blank())

png("false_positive_recrudescence_summary.png", height=2.5, width=8, units = "in", res=200)
show(fpr_summary_plot)
dev.off()

pdf("false_positive_recrudescence_summary.pdf", height=2.5, width=8)
show(fpr_summary_plot)
dev.off()



fpr_summary_plot_extended <- bind_rows(PfRecur_fpr %>% subset(epsilon==EPSILON_DEFAULT) %>% 
                                transmute(Site=Site, iter=iter, FPR=FPR, model=paste0("PfRecur\n(ω=", omega, ")")), 
                              match_counting_fpr, Plucinski_fpr) %>% 
  group_by(Site, model) %>% 
  summarise(median=median(FPR), 
            LCI=quantile(FPR, 0.025), UCI=quantile(FPR, 0.975),
            LCI2=quantile(FPR, 0.25), UCI2=quantile(FPR, 0.75)) %>% 
  #mutate(model=factor(model, levels=c("PfRecur", "CDC", "≥4/7", "≥5/7"))) %>% 
  ggplot(aes(x=model, y=median)) + 
  geom_point(cex=2) +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.55, col="#737373") + 
  geom_errorbar(aes(ymin=LCI2, ymax=UCI2), width=0.3, lwd=0.65) + 
  facet_grid(rows=vars(Site), scale="free_y") +
  scale_y_continuous(breaks=seq(0, 0.3, 0.02)) +
  xlab("Classification method") +
  ylab("False positive recrudescence rate") +
  theme_bw() + theme(panel.grid.major.x = element_blank())

png("false_positive_recrudescence_summary_extended.png", height=5, width=7, units = "in", res=200)
show(fpr_summary_plot_extended)
dev.off()

permuted_posterior_probs <- ggplot(PfRecur_permuted_data) + 
  geom_line(aes(x=M1, y=1-..y.., lwd="PfRecur", col=omega, group=omega), stat='ecdf') + 
  geom_line(aes(x=Plucinski, y=1-..y.., lwd="CDC model"), stat='ecdf', col="darkred") + 
  facet_grid(col=vars(Site)) + 
  xlab("Posterior probability of recrudescence") +
  ylab("Complementary ECDF") +
  coord_cartesian(ylim=c(0, 0.25)) + 
  scale_linewidth_manual(name="Model", values=c(0.9, 0.5)) + 
  scale_color_viridis_c(name="Per-clone\ndetection\nprob ω\n(PfRecur)") + 
  theme_bw() + theme(panel.spacing = unit(1.5, "lines"))

png("permuted_posterior_probs.png", height=3.2, width=10, units = "in", res=200)
show(permuted_posterior_probs)
dev.off()
