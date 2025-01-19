library(readr)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(reshape2)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(survival)
library(purrr)

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

# ============================ INPUT DATA ===============================

Angola_2021_data <- read_rds("../Angola_TES_data/Dimbu_Angola_TES_data.rds")
Angola_2021_pairs <-  read_rds("../Angola_TES_data/Dimbu_Angola_TES_pairs.rds")
Angola_2021_followup <- read_rds("../Angola_TES_data/Dimbu_Angola_TES_followup.rds")

MARKERS <- c(M313="M313", M383="M383", TA1="TA1", POLYA="POLYA",
             PFPK2="PFPK2", M2490="M2490", TA109="TA109") 

OMEGA_DEFAULT <- 0.9
EPSILON_DEFAULT <- 0.05
BETA <- 0.25

OMEGA_VALS <- c(0.75, 0.8, 0.85, seq(0.9, 1, 0.01))
EPSILON_VALS <- c(seq(0.01, 0.1, 0.01), 0.15, 0.2, 0.25)

MARKER_SETS <- list("7NMS"=MARKERS, "Exc_TA109"=setdiff(MARKERS, "TA109"))
MARKER_SET_NAMES <- c("7NMS"="All markers", "Exc_TA109"="Excluding TA109")
MARKER_SET_NAMES_long <- c("7NMS"="Classification\nbased on all 7 markers", 
                           "Exc_TA109"="Classification\nexcluding TA109")

RUN_ALL <- FALSE

Plucinski_classifications <- read_rds("../Plucinski_et_al_model/Full_dataset/Plucinski_classifications.rds")

Plucinski_params <- lapply(Plucinski_classifications, function(x) 
  x[["summary_statisics"]] %>% melt(id="Site") %>% 
    transmute(X=variable, site=Site, name=substr(site, 1, 1), value=value) %>% 
    extract(value, c("estimate", "LCI", "UCI"), "([0-9\\.]+) \\(([0-9\\.]+)–([0-9\\.]+)\\)")) %>% 
  bind_rows(.id="type")


# ============================ PARSE DATA ===============================

set.seed("20042024")

day_0_isolates <- 
  split(subset(Angola_2021_data$Sample_ID, Angola_2021_data$Timepoint=="D0"), 
        subset(substr(Angola_2021_data$Site, 1, 1), Angola_2021_data$Timepoint=="D0"))

Angola_2021_marker_set <- lapply(MARKERS, function(m) {
  sort(unique(setdiff(unlist(strsplit(Angola_2021_data[,m], "/")), "NA")))})

# rows: true allele, columns: observed allele
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

Angola_2021_MOI <- lapply(Angola_2021_genotype_matrix, function(x) {
  y <- Angola_2021_data[rownames(x), "MOI"]
  names(y) <- rownames(x)
  return(y)}) 

Angola_2021_locus_cardinality <- lapply(Angola_2021_genotype_matrix, rowSums)

# ============================ CLASSIFY RECURRENCES ============================

#PfRecur::evaluate_posterior(recurrent, ref_C, ref_I, genotype_matrix, error_matrix)

recurrent_list <- as.list(Angola_2021_pairs$D_REC)
ref_C_list <- as.list(Angola_2021_pairs$D_0)
ref_I_list <- lapply(ref_C_list, function(x) setdiff(day_0_isolates[[substr(x, 1, 1)]], x))

names(ref_C_list) <- Angola_2021_pairs$D_REC

# this iterates over 364 combinations of OMEGA, EPSILON and KEEP_MARKERS so can
# take ~15 minutes to run
if (!file.exists("Angola_TES_classifications.rds") | RUN_ALL) {
  Angola_2021_classifications <- lapply(EPSILON_VALS, function(EPSILON) {
      error_matrix <- PfRecur::geometric_error_matrix(Angola_2021_marker_set, EPSILON)
      prob_classifications <- lapply(names(MARKER_SETS), function(m) {
        parallel::mcmapply(PfRecur::evaluate_posterior, 
                             recurrent_list, ref_C_list, ref_I_list, 
                             MoreArgs=list(genotype_matrix=Angola_2021_genotype_matrix, 
                                           error_matrix=error_matrix, omega_vals = OMEGA_VALS, beta=BETA,
                                           keep_markers = MARKER_SETS[[m]]), SIMPLIFY = FALSE) %>%
          purrr::pmap(rbind) %>% lapply(function(x) x %>% as.data.frame %>% mutate(epsilon=EPSILON, markers=m))
        }) %>% purrr::pmap(rbind)
      print(EPSILON)
      return(prob_classifications)}) %>% purrr::pmap(rbind)
  write_rds(Angola_2021_classifications, "Angola_TES_classifications.rds")
} else {
  Angola_2021_classifications <- read_rds("Angola_TES_classifications.rds")
}

# ======================= VISUALISATION HELPER FUNCS ===========================
genotype_table <- function(my_code) {
  my_genotypes <- 
    Angola_2021_data[as.character(Angola_2021_pairs[my_code, c("D_0", "D_REC")]), MARKERS] %>%
    mutate_all(function(x) gsub("/", ", ", x )) %>% t #%>%
    #apply(1:2, function(x) paste0(strwrap(x, width=30), collapse="\n"))
  
  colnames(my_genotypes) <- paste0(c("Baseline", "Recurrent"), "\n(", 
                                   Angola_2021_pairs[my_code, c("D_0", "D_REC")], 
                                   ")\nParasitemia: ",
                                   formatC(c(Angola_2021_followup[my_code, "D0_parasitemia"],
                                             Angola_2021_followup[my_code, "DX_parasitemia"]),
                                           format = "e", digits=2),
                                   "/μL\nTemperature: ", 
                                   c(Angola_2021_followup[my_code, "D0_temperature"],
                                     Angola_2021_followup[my_code, "DX_temperature"]))
  return(my_genotypes)
}

genotype_table_basic <- function(my_code) {
  my_genotypes <- 
    Angola_2021_data[as.character(Angola_2021_pairs[my_code, c("D_0", "D_REC")]), MARKERS] %>%
    mutate_all(function(x) gsub("/", ", ", x )) %>% t #%>%
  #apply(1:2, function(x) paste0(strwrap(x, width=30), collapse="\n"))
  
  colnames(my_genotypes) <- paste0("Day ", c(0, Angola_2021_followup[my_code, "Followup"]))

    return(my_genotypes)
}

# ============================ SENSITIVITY TO FREE PARAMS ============================

free_param_sens <- Angola_2021_classifications[["metrics"]] %>% 
  split(f=.$recurrence) %>% lapply(function(x) {
    
    recurrence <- as.character(x[1, "recurrence"])
    paired_baseline <- ref_C_list[[recurrence]]
    study_arm <- substr(recurrence, 1, 1)
    
    q_est <- Plucinski_params %>% subset(name==study_arm & X=="q")
    d_est <- Plucinski_params %>% subset(name==study_arm & X=="d")
    
    prob_recur_plot <- split(x, f=x$markers) %>% lapply(function(y) {
      y <- as.data.frame(y)
      q_est <- Plucinski_params %>% subset(name==study_arm & X=="q" & type==y[1, "markers"])
      d_est <- Plucinski_params %>% subset(name==study_arm & X=="d" & type==y[1, "markers"])
      p <- y %>% select(-recurrence, -markers) %>% reshape2::melt(id=c("omega", "epsilon")) %>%
        ggplot(aes(x=as.character(omega), y=as.character(epsilon))) + 
        geom_tile(aes(fill=round(value, 5)), col="black") + 
        geom_text(aes(label=round(value, 2), 
                      size=(omega==OMEGA_DEFAULT | epsilon==EPSILON_DEFAULT),
                      fontface=ifelse(omega==OMEGA_DEFAULT | epsilon==EPSILON_DEFAULT, 
                                      "bold", "plain"))) + 
        facet_grid(cols=vars(variable)) + 
        scale_fill_gradient(name="Posterior\nmetric)", 
                            low="#f2c4a2", high="darkred", limits=c(0, 1)) + 
        scale_size_manual(values=c(1.7, 2.05)) +
        guides(size="none") +
        xlab(bquote(omega ~ "(per-clone marker-wise probability of detection)")) + 
        ylab(bquote(epsilon ~ "(genotyping error rate)")) +
        ggtitle(MARKER_SET_NAMES[y[1, "markers"]]) +
        labs(caption=paste0("Posterior probability of recrudescence based on CDC model: ",
                            round(Angola_2021_pairs[substr(recurrence, 1, 8), y[1, "markers"]], 3), "\n",
                            "Estimated per-clone unobserved allele penalty: ", q_est[,"estimate"], " (", q_est[,"LCI"], "-", q_est[,"UCI"], ")\n",
                            "Estimated genotyping error rate: ", round(1-as.numeric(d_est[,"estimate"]), 3), " (", 
                            round(1-as.numeric(d_est[,"UCI"]), 3), "-", round(1-as.numeric(d_est[,"LCI"]), 3), ")")) + 
        theme_bw() + theme(axis.text.x.bottom = element_text(size=7),
                           axis.title.x.top = element_text(face="italic"), 
                           axis.title.y.right = element_text(face="italic"), 
                           axis.text.y = element_text(size=7), 
                           plot.title = element_text(hjust=0.5, face="bold", size=11.5),
                           plot.subtitle = element_text(hjust=0.5, face="italic"),
                           plot.caption = element_text(hjust=0.5, color="#545454", size=10),
                           legend.position = "none")
      return(p)})
    
    genotype_summary <- tableGrob(genotype_table(substr(recurrence, 1, 8)), theme=ttheme_default(base_size = 10.5))
    genotype_summary$widths <- unit(c(1/5, 2/5, 2/5), "npc")
  
    return(plot_grid(genotype_summary, plot_grid(plotlist=prob_recur_plot), nrow=1, 
                    rel_widths = c(0.8, 2), scale=0.925, align="vh", axis="tblr"))
})

if (!dir.exists("Free_param_senstivity/")) dir.create("Free_param_sensitivity/")

lapply(names(free_param_sens), function(x) {
  png(paste0("Free_param_sensitivity/", x, ".png"), width=20, height=4, unit="in", res=125)
  show(free_param_sens[[x]])
  dev.off()
})


# ============================ SENSITIVE TO OMEGA ============================

omega_sens <- c("LL21-061", "ZL21-233", "LL21-010", "ZQ21-085", "LL21-077", "BD21-053", # stronger dependence  
                "ZQ21-103", "ZL21-245", #flips when zero vs nonzero
                 "ZL21-203", "ZL21-269", "ZL21-304", "LL21-054", "ZL21-260") #clear recrudescence but diff mixtures

omega_sens_plots <- lapply(omega_sens, function(i) {

  prob_recur_plot <- Angola_2021_classifications[["metrics"]] %>% 
    subset(recurrence==Angola_2021_pairs[i, "D_REC"] & epsilon==EPSILON_DEFAULT) %>%
    split(f=.$markers) %>% lapply(function(x) {
      x <- as.data.frame(x)
      markers <- x[1, "markers"]
      q_est <- Plucinski_params %>% subset(name==substring(Angola_2021_pairs[i, "D_0"], 1, 1) & type==markers & X=="q")
      d_est <- Plucinski_params %>% subset(name==substring(Angola_2021_pairs[i, "D_0"], 1, 1) & type==markers & X=="d")
      p <- x %>% select(-recurrence, -epsilon, -markers) %>%
        reshape2::melt(id="omega") %>%
        mutate(Metric=ifelse(variable=="M1", "M1 (at least one recrudescent clone)", 
                             "M2 (proportion of recrudescent clones)")) %>%
        ggplot(aes(x=omega, y=value, color=Metric, pch=Metric, group=Metric, lty=Metric)) + 
        geom_hline(aes(yintercept=round(Angola_2021_pairs[i, x[1, "markers"]], 3)),
                   col="#545454", lwd=1, alpha=0.5) +
        geom_line(lwd=0.4) + geom_point(cex=2.2) + 
        xlab(bquote(omega ~ "(per-clone marker-wise probability of detection)")) +
        scale_y_continuous(limits=c(-0.01, 1.01), name = "Posterior estimate") +
        scale_linetype_manual(values=c(1, 1)) + 
        scale_color_manual(values=c("navyblue", "#6dc9ba")) +
        scale_shape_manual(values=c(16, 13)) +
        labs(subtitle=paste0("Genotyping error rate (normalised geometric model): ", EPSILON_DEFAULT)) +
        ggtitle(MARKER_SET_NAMES[markers]) +
        labs(caption=paste0("Posterior probability of recrudescence based on CDC model: ",
                            round(Angola_2021_pairs[i, markers], 3), "\n",
                            "Estimated per-clone unobserved allele penalty: ", q_est[,"estimate"], " (", q_est[,"LCI"], "-", q_est[,"UCI"], ")\n",
                            "Estimated genotyping error rate: ", round(1-as.numeric(d_est[,"estimate"]), 3), " (", 
                            round(1-as.numeric(d_est[,"UCI"]), 3), "-", round(1-as.numeric(d_est[,"LCI"]), 3), ")")) + 
        theme_bw() + theme(plot.title = element_text(hjust=0.5, face="bold", size=11.5),
                           plot.subtitle = element_text(hjust=0.5, face="italic"),
                           plot.caption = element_text(hjust=0.5, color="#545454", size=10),
                           legend.text = element_text(size=8),
                           legend.title = element_blank(),
                           legend.box.spacing = unit(1, "pt"),
                           legend.position = "top")
      return(p)})
  
  genotype_summary <- tableGrob(genotype_table(i), theme=ttheme_default(base_size = 10.5))
  genotype_summary$widths <- unit(c(1/5, 2/5, 2/5), "npc")
  
  return(plot_grid(genotype_summary, plot_grid(plotlist=prob_recur_plot), nrow=1, 
                   rel_widths = c(1, 2), scale=0.925, align="vh", axis="tblr"))
})
names(omega_sens_plots) <- omega_sens


my_title_1 <- title <- ggdraw() + 
  draw_label("(A) Classification sensitive to per-clone marker-wise probability of detection",
             fontface = 'bold.italic', x = 0.5, hjust = 0.5, size=19) 
my_title_2 <- title <- ggdraw() + 
  draw_label("(B) Accommodating imperfect detection protects against possible human error",
             fontface = 'bold.italic', x = 0.5, hjust = 0.5, size=19) 
my_title_3 <- title <- ggdraw() + 
  draw_label("(C) Resolution of mixtures sensitive to per-clone marker-wise probability of detection",
             fontface = 'bold.italic', x = 0.5, hjust = 0.5, size=19) 

my_plots_1 <- plot_grid(my_title_1, plot_grid(plotlist=omega_sens_plots[1:6], 
                                              ncol=1, labels=paste0("(A", 1:2, ")")), 
                        ncol=1, rel_heights = c(0.025, 1.2))

my_plots_2 <- plot_grid(my_title_2, plot_grid(plotlist=omega_sens_plots[7:8], 
                                              ncol=1, labels=paste0("(B", 1:2, ")")), 
                        ncol=1, rel_heights = c(0.025, 0.4))

my_plots_3 <- plot_grid(my_title_3, plot_grid(plotlist=omega_sens_plots[9:12], 
                                              ncol=1, labels=paste0("(C", 1:4, ")")), 
                        ncol=1, rel_heights = c(0.025, 0.8))

png("Dimbu_Angola_latent_baseline_1.png", width=16, height=21, units="in", res=180)
show(my_plots_1)
dev.off()

png("Dimbu_Angola_latent_baseline_2.png", width=16, height=21, units="in", res=180)
show(plot_grid(my_plots_2, my_plots_3, ncol=1, rel_heights = c(1, 2)))
dev.off()


# ============================ SENSITIVE TO EPSILON ============================

epsilon_sens <- c("ZL21-260", "ZQ21-085", "ZL21-233", "LL21-061", "BD21-053")

epsilon_sens_plots <- lapply(epsilon_sens, function(i) {
  
  prob_recur_plot <- Angola_2021_classifications[["metrics"]] %>% 
    subset(recurrence==Angola_2021_pairs[i, "D_REC"] & omega==OMEGA_DEFAULT) %>%
    split(f=.$markers) %>% lapply(function(x) {
      x <- as.data.frame(x)
      markers <- x[1, "markers"]
      q_est <- Plucinski_params %>% subset(name==substring(Angola_2021_pairs[i, "D_0"], 1, 1) & type==markers & X=="q")
      d_est <- Plucinski_params %>% subset(name==substring(Angola_2021_pairs[i, "D_0"], 1, 1) & type==markers & X=="d")
      p <- x %>% select(-recurrence, -omega, -markers) %>%
        reshape2::melt(id="epsilon") %>%
        mutate(Metric=ifelse(variable=="M1", "M1 (at least one recrudescent clone)", 
                             "M2 (proportion of recrudescent clones)")) %>%
        ggplot(aes(x=epsilon, y=value, color=Metric, pch=Metric, group=Metric, lty=Metric)) + 
        geom_hline(aes(yintercept=round(Angola_2021_pairs[i, x[1, "markers"]], 3)),
                   col="#545454", lwd=1, alpha=0.5) +
        geom_line(lwd=0.4) + geom_point(cex=2.2) + 
        xlab(bquote(epsilon ~ "(genotyping error rate)")) +
        scale_y_continuous(limits=c(-0.01, 1.01), name = "Posterior estimate") +
        scale_linetype_manual(values=c(1, 1)) + 
        scale_color_manual(values=c("navyblue", "#6dc9ba")) +
        scale_shape_manual(values=c(16, 13)) +
        labs(subtitle=paste0("Per-clone marker-wise probability of detection: ", OMEGA_DEFAULT)) +
        ggtitle(MARKER_SET_NAMES[markers]) +
        labs(caption=paste0("Posterior probability of recrudescence based on CDC model: ",
                            round(Angola_2021_pairs[i, markers], 3), "\n",
                            "Estimated per-clone unobserved allele penalty: ", q_est[,"estimate"], " (", q_est[,"LCI"], "-", q_est[,"UCI"], ")\n",
                            "Estimated genotyping error rate: ", round(1-as.numeric(d_est[,"estimate"]), 3), " (", 
                            round(1-as.numeric(d_est[,"UCI"]), 3), "-", round(1-as.numeric(d_est[,"LCI"]), 3), ")")) + 
        theme_bw() + theme(plot.title = element_text(hjust=0.5, face="bold", size=11.5),
                           plot.subtitle = element_text(hjust=0.5, face="italic"),
                           plot.caption = element_text(hjust=0.5, color="#545454", size=10),
                           legend.text = element_text(size=8),
                           legend.title = element_blank(),
                           legend.box.spacing = unit(1, "pt"),
                           legend.position = "top")
      return(p)})
  
  genotype_summary <- tableGrob(genotype_table(i), theme=ttheme_default(base_size = 10.5))
  genotype_summary$widths <- unit(c(1/5, 2/5, 2/5), "npc")
  
  return(plot_grid(genotype_summary, plot_grid(plotlist=prob_recur_plot), nrow=1, 
                   rel_widths = c(1, 2), scale=0.925, align="vh", axis="tblr"))
})
names(epsilon_sens_plots) <- epsilon_sens


my_title <- title <- ggdraw() + 
  draw_label("Classification sensitive to genotyping error",
             fontface = 'bold.italic', x = 0.5, hjust = 0.5, size=19) 

my_plots <- plot_grid(my_title, plot_grid(plotlist=epsilon_sens_plots, 
                                            ncol=1, labels=paste0("(A", 1:5, ")")), 
                        ncol=1, rel_heights = c(0.025, 1))

png("Dimbu_Angola_genotyping_error.png", width=16, height=17.5, units="in", res=180)
show(my_plots)
dev.off()

# ============================ MODEL COMPARISON ===============================

posterior_summary <- Angola_2021_classifications[["metrics"]] %>% 
  subset(epsilon==EPSILON_DEFAULT & omega==OMEGA_DEFAULT) %>% 
  transmute(pair=substr(recurrence, 1, 8), markers=markers, prob=M1, model="PfRecur")

Dimbu_summary <- Angola_2021_pairs %>% select(-D_0, -D_REC, -Dimbu_posterior) %>%
  reshape2::melt(id="Code") %>% 
  rename(pair=Code, markers=variable, prob=value) %>% 
  mutate(model="CDC")

isolate_ordering <- (Angola_2021_pairs %>% arrange(`7NMS`))$Code
posterior_summary <- bind_rows(posterior_summary, Dimbu_summary) %>% 
  merge(Angola_2021_followup) %>% subset(D0_MOI>1) %>%
  mutate(pair=factor(pair, levels=isolate_ordering),
         markers=MARKER_SET_NAMES_long[markers])

prob_plot <- reshape2::dcast(posterior_summary, pair+markers~model, value.var="prob") %>% 
  mutate(pair=factor(pair, levels=isolate_ordering)) %>%
  ggplot(aes(y=pair, x=PfRecur)) + 
  geom_point(aes(x=CDC, color=CDC, pch="P"), cex=3.5) + 
  geom_point(aes(color=PfRecur, pch="C"), cex=2.2) + 
  geom_segment(aes(xend=CDC, lty=CDC-PfRecur>0), 
               col="#575757", lwd=0.35) + 
  scale_linetype_manual(values=c(1, 2), guide="none") +
  scale_shape_manual(values=c(16, 21), name="Model", 
                     labels=c(paste0("Probability of at least one recrudescent clone (M1) under PfRecur\n(genotyping error rate: ", 
                                     EPSILON_DEFAULT, ", per-clone marker-wise detection prob: ", OMEGA_DEFAULT, ")"), 
                              "Posterior probability of recrudescence under CDC model")) +
  scale_color_gradient(low="#cc8d3b", high="darkred", guide="none") + 
  facet_grid(cols=vars(markers)) + 
  guides(shape = guide_legend(ncol = 1)) +
  xlab("Posterior probability of recrudescence") +
  theme_bw() + theme(axis.text.y = element_blank(),
                     legend.position = "bottom",
                     legend.title.position = "top", 
                     legend.title = element_text(hjust=0.5, face="italic"),
                     axis.ticks.length.y = unit(6, "pt"),
                     #panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank(),
                     #legend.box.spacing = unit(5, "pt"),
                     axis.title.y = element_blank())

metric_plot <- posterior_summary %>% 
  select(pair, D0_MOI, DX_MOI, D0_parasitemia, DX_parasitemia) %>% melt(id="pair") %>% 
  extract(variable, c("Timepoint", "Metric"), "(D[0X])_([A-Za-z]+)") %>% unique %>% 
  split(f=.$Timepoint) %>% lapply(function(x) 
    x %>% select(-Timepoint) %>% dcast(pair~Metric, value.var = "value")) %>% 
  bind_rows(.id="Timepoint") %>% mutate(pair=factor(pair, levels=isolate_ordering)) %>% 
  ggplot(aes(y=pair, x=1, label=MOI)) + 
  geom_tile(aes(fill=parasitemia), col="black") + geom_text(size=3) + 
  facet_grid(cols=vars(paste0(Timepoint, "\n MOI"))) + 
  scale_fill_gradient(trans="log10", low="white", high="darkred", name="Parasitemia (per μL)") +
  scale_x_continuous(expand=c(0, 0)) + 
  xlab("") + ylab("Patient ID") +
  theme_bw() + theme(axis.text.x=element_blank(), 
                     axis.ticks.x = element_blank(), 
                     #legend.box.spacing = unit(5, "pt"),
                     legend.title.position = "top", 
                     legend.text = element_text(angle=90),
                     legend.title = element_text(hjust=0.5, face="italic"),
                     legend.position="bottom") 

summary_plot <- plot_grid(metric_plot, prob_plot, align="h", axis="tblr", rel_widths = c(1, 2), nrow=1)

png("Dimbu_Angola_summary.png", width=7, height=8, units="in", res=200)
show(summary_plot)
dev.off()

# ============================ DIFFERENT PREDICTIONS ============================

highlight_isolates <- 
  unique(as.character((dcast(posterior_summary, pair+markers~model, value.var="prob") %>% 
                         subset(abs(PfRecur-CDC)>0.1) %>%
                         arrange(-abs(PfRecur-CDC)))$pair))

metric_names <- c("M1"=">=1 clone\nrecrudescent",
                  "M2"="Proportion\nrecrudescent",
                  "Plucinski"="Probability\nrecrudescence")

MARKER_SET_NAMES <- c("7NMS"="All markers", "Exc_TA109"="Exc. TA109")
prob_palette <- colorRampPalette(c("#cc8d3b", "darkred"))(1001)

model_diff_plots <- lapply(highlight_isolates, function(i) {

  PfRecur_prob <- Angola_2021_classifications[["metrics"]] %>% as.data.frame %>%
    subset(recurrence==Angola_2021_pairs[i, "D_REC"] & omega==OMEGA_DEFAULT & epsilon==EPSILON_DEFAULT) %>%
    select(-omega, -recurrence, -epsilon) %>%
    reshape2::melt(id="markers") %>% mutate(model="PfRecur") %>% mutate(pair=i)
  
  Plucinski_prob <- Dimbu_summary %>% subset(pair==i) %>% 
    rename(value=prob) %>% mutate(variable="Plucinski")
  
  summary_prob <- bind_rows(PfRecur_prob, Plucinski_prob) %>% 
    mutate(value=round(value, 3), markers=MARKER_SET_NAMES[markers]) %>% 
    dcast(variable~markers, value.var="value") %>% 
    mutate(variable=metric_names[variable]) %>% 
    tibble::column_to_rownames(var="variable") %>% t 
  
  prob_cols <- sapply(summary_prob, function(x) prob_palette[round(x*1000)+1])
  
  summary_prob <- tableGrob(summary_prob, 
                            theme=ttheme_default(base_size=10.5,
                                                 colhead=list(bg_params=list(fill=NA),
                                                              fg_params=list(col="#545454", fontface=1L)),
                                                 rowhead=list(fg_params=list(col="#545454")),
                                                 core=list(bg_params=list(fill=prob_cols, alpha=0.6))))

  model_header <- tableGrob(matrix(0, nrow=1, ncol=3), rows="", 
                            cols=c("PfRecur M1", "PfRecur M2", "CDC"),
                            theme=ttheme_minimal(base_size = 10.5,
                                                 colhead=list(fg_params=list(col="#545454"))))
  
  prob_table <- gtable_combine(model_header[1,], summary_prob, along=2)
  
  # https://stackoverflow.com/questions/15622001/how-to-display-only-integer-values-on-an-axis-using-ggplot2
  genotype_summary <- tableGrob(genotype_table_basic(i), theme=ttheme_default(base_size = 10.5))
  genotype_summary$widths <- unit(c(1/6, 5/12, 5/12), "npc")
  
  pair_title <- ggdraw() + draw_label(i, fontface = 'bold', x = 0.5, hjust = 0.5, size=12)
  
  return(plot_grid(pair_title, prob_table, genotype_summary, rel_heights = c(0.35, 1, 2.2), ncol=1))
})
names(model_diff_plots) <- highlight_isolates

# ============================ IMPLICATIONS =================================

# code adapted from https://github.com/MateuszPlucinski/AngolaTES2021
KM_estimate <- function(followup, outcome, duration) {
  outcome_summary <- Surv(followup, outcome)
  KM_fit <- survfit(outcome_summary~1, conf.lower="peto")
  KM_summary <- summary(KM_fit, times=duration)[c("surv", "lower", "upper")]
  return(unlist(KM_summary))
}

# code adapted from https://github.com/MateuszPlucinski/AngolaTES2021
KM_probabilistic <- function(arm_summary, prob_colname, niter) {
  n_indiv <- nrow(arm_summary)
  prob_recr <- arm_summary[,prob_colname]
  followup <- arm_summary$Followup
  duration <- max(followup)
  
  if (max(prob_recr)>0) {
    KM_summary <- sapply(1:niter, function(i) 
      KM_estimate(followup, as.numeric(runif(n_indiv)<prob_recr), duration))
    return(rowMeans(KM_summary))
  } else {
    return(c(surv=1, lower=1, upper=1))
  }
}


NITER <- 10000

TES_summary <- Angola_2021_classifications[["metrics"]] %>% as.data.frame %>%
  subset(epsilon==EPSILON_DEFAULT & omega==OMEGA_DEFAULT & markers=="7NMS") %>% 
  mutate(pair=substr(recurrence, 1, 8)) %>% 
  select(-recurrence, -omega, -epsilon, -markers) %>%
  merge(Angola_2021_followup, all=TRUE) %>%
  merge(Plucinski_classifications[["7NMS"]][["posterior_recrudesence"]] %>% transmute(pair=isolate, Plucinski_7NMS=prob), all=TRUE) %>%
  mutate(arm=paste0(Province, "_", Drug))
TES_summary[is.na(TES_summary)] <- 0

TES_summary_by_arm <- split(TES_summary, TES_summary$arm)

efficacy_summary <- lapply(TES_summary_by_arm, function(arm_summary) 
  sapply(c("M1", "Dimbu_posterior", "Plucinski_7NMS"), function(i) KM_probabilistic(arm_summary, i, NITER)))

write_rds(efficacy_summary, "Angola_TES_efficacy_summary.rds")

png("Dimbu_Angola_diff.png", width=28, height=12, units="in", res=180)
plot_grid(plotlist=model_diff_plots, nrow=3, 
          labels=paste0("(", LETTERS[1:length(highlight_isolates)], ")"))
dev.off()
