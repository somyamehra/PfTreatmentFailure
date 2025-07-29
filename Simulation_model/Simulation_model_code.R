library(tidyr)
library(readr)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(ggridges)
library(parallel)

library(showtext)
font_add(family = "Arial", regular = "Arial.ttf") ## here is the path to the font to add.
showtext.auto()

library(extrafont)
font_import() # takes a few minutes
loadfonts(device="postscript")

# https://rstats-tips.net/2020/07/31/get-rid-of-info-of-dplyr-when-grouping-summarise-regrouping-output-by-species-override-with-groups-argument/
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

library(PfRecur)
library(PDQutils)
library(copula)
library(poisbinom)
library(VGAM)

setwd(this.path::here())

# ============================ SIMULATION MODEL ===============================

simulate_data <- function(n_marker, marker_cardinality, partition_prob, 
                          error_prob, mean_MOI, max_MOI, omega,
                          n_baseline, n_datasets) {
  # linear shape
  pop_allele_freq <- lapply(marker_cardinality, function(x) (1:x)/(x*(x+1)/2))
  
  # structured genotyping error (allele frequencies reordered for each dataset)
  delta <- lapply(marker_cardinality, function(m) {
    error_matrix <- matrix(1, nrow=m, ncol=m, dimnames = list(1:m, 1:m))
    if (m==1) return(error_matrix)
    diag(error_matrix) <- 1 - error_prob
    if (error_prob>0) {
      for (i in 1:m) {
        non_self <- setdiff(1:m, i)
        p <- dgeom(abs(non_self-i), 1-error_prob)
        error_matrix[i,non_self] <- p*error_prob/sum(p)
      }
    }
    return(error_matrix)
  })
  names(delta) <- 1:n_marker
  
  observed_allele <- function(true_allele, error_matrix) {
    sapply(true_allele, function(x) sample(1:nrow(error_matrix), 1, prob=error_matrix[x,]))
  }
  
  sibling_partition <- function(MOI) {
    table(cumsum(sample(c(0, 1), size=MOI, replace=TRUE, 
                        prob=c(partition_prob, 1-partition_prob))))
  }
  
  sample_siblings <- function(n, pop_allele_freq_reordered) {
    parents <- sapply(pop_allele_freq_reordered,
                      function(x) sample(1:length(x), 2, replace=TRUE, prob=x))
    siblings <- apply(parents, 2, function(x) sample(x, n, replace=TRUE))
    if (n==1) return(t(as.matrix(siblings)))
    return(siblings)
  }
  
  sample_isolate <- function(MOI, pop_allele_freq_reordered) {
    do.call(rbind, sapply(sibling_partition(MOI), sample_siblings, pop_allele_freq_reordered, simplify = FALSE))
  }
  
  generate_recurrence <- function(paired_clones, n_recrudescent, final_MOI, pop_allele_freq_reordered) {
    if (n_recrudescent>=1) {
      recrudescent_clones <- 
        paired_clones[sample(1:nrow(paired_clones), n_recrudescent, replace=FALSE),,drop=FALSE]
    }
    
    if (final_MOI-n_recrudescent>=1) {
      new_clones <- sample_isolate(final_MOI-n_recrudescent, pop_allele_freq_reordered)
    }
    
    if (n_recrudescent==0) {
      recurrent_clones <- new_clones
    } else if (final_MOI-n_recrudescent==0) {
      recurrent_clones <- recrudescent_clones
    } else {
      recurrent_clones <- rbind(recrudescent_clones, new_clones)
    }
    
    recurrent_obs_allele <- 
      apply(recurrent_clones, 2, function(x) unique(x[rbinom(length(x), 1, omega)>0]),
            simplify = FALSE)
    #recurrent_obs_allele <- apply(recurrent_clones, 2, unique, simplify = FALSE)
    #mapply(observed_allele, apply(recurrent_clones, 2, unique), delta)
    
    return(recurrent_obs_allele)
  }
  
  # MOI sampled from truncated geometric [1, MAX_MOI]
  MOI_dist <- dgeom(0:(max_MOI-1), 1/mean_MOI)/pgeom(max_MOI-1, 1/mean_MOI)
  
  simulated_data <- lapply(1:n_datasets, function(k) {
    print(k)
    
    pop_allele_freq_reordered <- lapply(pop_allele_freq, function(x) sample(x, length(x)))
    
    baseline_MOI <- sample(1:max_MOI, n_baseline, replace=TRUE, prob=MOI_dist)

    # column = marker, row = clone
    baseline_clones <- lapply(baseline_MOI, sample_isolate, pop_allele_freq_reordered)
    
    # adjusted for error; row = isolate, column = allele
    baseline_obs_genotype_matrix <- 
      lapply(marker_cardinality, function(x) matrix(0, nrow=n_baseline, ncol=x))
    
    for (isolate in 1:n_baseline) {
      for (marker in 1:n_marker) {
        detected_clones <- rbinom(nrow(baseline_clones[[isolate]]), 1, omega)>0
        if (any(detected_clones)) {
          obs_alleles <- observed_allele(unique(baseline_clones[[isolate]][,marker][detected_clones]), 
                                         delta[[marker]])
          baseline_obs_genotype_matrix[[marker]][isolate, obs_alleles] <- 1
        }
      }
    }
    
    # assumed MOI = maximum cardinality at any locus
    baseline_obs_MOI <- apply(sapply(baseline_obs_genotype_matrix, rowSums), 1, max)
    
    genotype_matrix <- lapply(baseline_obs_genotype_matrix, function(x) {
      x <- rbind(x, 0)
      rownames(x) <- c(paste0("B", 1:n_baseline), "R")
      colnames(x) <- 1:ncol(x)
      return(x) })
    names(genotype_matrix) <- 1:n_marker

    
    recurrence_summary <- function(n_recrudescent, recurrent_MOI, pair) {
      
      target_alleles <- generate_recurrence(baseline_clones[[pair]], n_recrudescent, 
                                            recurrent_MOI, pop_allele_freq_reordered)
      
      recurrent_obs_MOI <- max(sapply(target_alleles, length))
      
      if (n_recrudescent>recurrent_obs_MOI) return(NULL)
      
      for (i in 1:n_marker) {
        genotype_matrix[[i]]["R", ] <- 0
        genotype_matrix[[i]]["R", as.character(target_alleles[[i]])] <- 1
      }
      
      genotype_matrix <- lapply(genotype_matrix, function(x) x[rowSums(x)>0,])
      
      posterior <- 
        PfRecur::evaluate_posterior("R", paste0("B", pair), paste0("B", setdiff(1:n_baseline, pair)),
                                    genotype_matrix, delta, omega, beta=1)
      

      return(c(pair=pair, n_recrudescent=n_recrudescent, 
               recurrent_MOI=recurrent_MOI,
               baseline_MOI=baseline_MOI[pair],
               recurrent_obs_MOI=recurrent_obs_MOI,
               baseline_obs_MOI=baseline_obs_MOI[pair],
               M1=as.numeric(posterior[["metrics"]][1, "M1"]),
               M2=as.numeric(posterior[["metrics"]][1, "M2"])))
    }
    
    classification_summary <- mclapply(which(baseline_obs_MOI>0), function(pair) {
      lapply(1:max_MOI, function(x) {
        lapply(0:min(x, baseline_obs_MOI[pair]), recurrence_summary, x, pair) %>%
          Filter(function(x) !is.null(x), .) %>% do.call(rbind, .)}) %>% 
        do.call(rbind, .) %>% as.data.frame}) %>% bind_rows()
    
    return(classification_summary)}) %>% bind_rows(.id="iter")
  
  return(list(simulated_data=simulated_data,
              omega=omega,
              n_marker=n_marker, marker_cardinality=marker_cardinality,
              error_prob=error_prob, mean_MOI=mean_MOI, max_MOI=max_MOI,
              n_baseline=n_baseline, partition_prob=partition_prob))
}

if (!file.exists("Simulation_results.rds")) {
  set.seed("6216")
  sim_data <- simulate_data(n_marker=7, marker_cardinality=c(30, 30, 20, 20, 20, 10, 10), 
                            partition_prob=0.1, error_prob=0.05, mean_MOI=3, 
                            max_MOI=9, omega=0.9, n_baseline=25, n_dataset=40)
  write_rds(sim_data, "Simulation_results.rds", compress = "gz")
} else {
  sim_data <- read_rds("Simulation_results.rds")
}


conservative_plot <- ggplot(sim_data[[1]] %>% subset(recurrent_MOI%%2==1 & n_recrudescent<9) %>% 
         mutate(recurrent_MOI_label=paste0("Recurrent\nMOI: ", recurrent_MOI))) + 
  geom_density_ridges(aes(x=M2, y=paste0(n_recrudescent, "/", recurrent_MOI)), 
                      stat="binline", bins=50, size=0.25, alpha=0.25, fill="darkorange") +
  geom_segment(aes(x=n_recrudescent/recurrent_MOI, xend=n_recrudescent/recurrent_MOI, 
                   y=n_recrudescent+0.95, yend=n_recrudescent+2), lwd=0.6, col="brown") + 
  facet_grid(rows=vars(recurrent_MOI_label), scale="free_y", space="free_y") + 
  xlab("Posterior proportion of recrudescent clones") +
  ylab("True proportion of recrudescent clones") +
  ggtitle("M2: proportion of recrudescent clones") +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(hjust=0.5))

anticonservative_metric <- sim_data[[1]] %>% 
  mutate(type=factor(ifelse(n_recrudescent==0, "0 recrudescent clones", "\u2265 1 recrudescent clone"), 
                     levels=c("0 recrudescent clones", "\u2265 1 recrudescent clone")))

sens_spec <- anticonservative_metric %>% group_by(recurrent_MOI, type) %>% 
  summarise(metric=mean(M1>0.5)) %>%
  mutate(label=ifelse(!type=="0 recrudescent clones", paste0("Sensitivity: ", round(metric, 3)), 
                      paste0("Specificity: ", round(1-metric, 3)))) %>%
  mutate(type=factor(type, levels=c("0 recrudescent clones", "\u2265 1 recrudescent clone")))

anticonservative_plot <- ggplot() + 
  geom_density_ridges(data=anticonservative_metric,
                      aes(x=M1, y=paste0(recurrent_MOI), fill=type),
                      stat="binline", bins=30, size=0.25, alpha=0.25, fill="darkorange") +
  geom_segment(data=anticonservative_metric,
               aes(x=as.numeric(n_recrudescent!=0), xend=as.numeric(n_recrudescent!=0), 
                   y=recurrent_MOI-0.05, yend=recurrent_MOI+0.8), lwd=0.6, col="brown") +
  annotate("segment", x=0.5, xend=0.5, y=-Inf, yend=Inf, col="navyblue", lty=2) +
  geom_label(data=sens_spec, aes(x=0.5, y=recurrent_MOI, label=label), vjust=-0.6, size=3) +
  facet_grid(cols=vars(type)) +
  xlab("Posterior probability \u2265 1 recrudescent clone") +
  ylab("Recurrent MOI") +
  ggtitle("M1: \u2265 1 recrudescent clone") +
  scale_x_continuous(breaks=seq(0, 1, 0.25)) +
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(hjust=0.5),
                     text=element_text(family="sans"))

png("Simulation_results.png", height=9, width=9, units = "in", res=200)
plot_grid(anticonservative_plot, NULL, conservative_plot, nrow=1, 
          rel_widths = c(1, 0.025, 1), labels=c("(A)", "", "(B)"))
dev.off()

pdf("Simulation_results.pdf", height=9, width=9)
plot_grid(anticonservative_plot, NULL, conservative_plot, nrow=1, 
          rel_widths = c(1, 0.025, 1), labels=c("(A)", "", "(B)"))
dev.off()



