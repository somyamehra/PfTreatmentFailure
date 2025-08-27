library(readxl)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(reshape2)
library(tidyr)
library(stringr)
library(copula)
library(ggplot2)
library(cowplot)
library(readr)

library(showtext)
font_add(family = "Arial", regular = "Arial.ttf") ## here is the path to the font to add.
showtext.auto()

library(extrafont)
font_import() # takes a few minutes
loadfonts(device="postscript")

setwd(this.path::here())

# ============================ HELPER FUNCTIONS ============================

add_prefix <- function(allele, prefix) {
  ifelse(!is.na(allele), gsub(", ", paste0("/", prefix), paste0(prefix, allele)), NA)
}

calc_MOI <- function(allele) {
  ifelse(allele!="", str_count(allele, "/")+1, 0)
}

expected_prop_IBS <- function(MOI, card) {
  if (MOI*card==0) return(NaN)
  return(sum(sapply(1:(MOI-card+1), function(k) choose(MOI, k)*Stirling2(k, 1)*
                      Stirling2(MOI-k, card-1)/Stirling2(MOI, card)) * 
               ((1:(MOI-card+1))/MOI)^2))
}


# ============================ PARSE GENOTYPES ============================

Angola_2021 <- bind_rows(read_xlsx("Dimbu_Angola_TES_genotypes.xlsx", skip=3),
                         read_xlsx("Dimbu_Angola_TES_genotypes.xlsx", sheet=2, skip=3)) %>%
  as.data.frame()

MARKERS <- c(M313="M313", M383="M383", TA1="TA1", POLYA="POLYA",
             PFPK2="PFPK2", M2490="M2490", TA109="TA109") 

TYPES <- c("7NMS", "Exc_TA109")
STUDY_SITES <- c("Benguela", "Lunda Sul", "Zaire")

Angola_2021_data <- Angola_2021 %>% rename(Sample_ID=`Sample.ID`) %>%
  extract(Sample_ID, c("Code", "Timepoint"), "([A-Z]{2}[0-9-]{6,})(D[0-9]+)", remove = F) %>%
  unite(TA109, colnames(Angola_2021)[grepl("TA109", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  unite(M313, colnames(Angola_2021)[grepl("313", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  unite(M383, colnames(Angola_2021)[grepl("383", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  unite(TA1, colnames(Angola_2021)[grepl("TA1_", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  unite(POLYA, colnames(Angola_2021)[grepl("POLYA", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  unite(PFPK2, colnames(Angola_2021)[grepl("PFPK2", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  unite(M2490, colnames(Angola_2021)[grepl("2490", colnames(Angola_2021))], sep=", ", na.rm=T) %>%
  mutate(TA109=add_prefix(TA109, ""),
         M313=add_prefix(M313, ""),
         M383=add_prefix(M383, ""),
         TA1=add_prefix(TA1, ""),
         POLYA=add_prefix(POLYA, ""),
         PFPK2=add_prefix(PFPK2, ""),
         M2490=add_prefix(M2490, "")) %>%
  mutate(MOI=pmax(calc_MOI(TA109), calc_MOI(M313), calc_MOI(M383), calc_MOI(TA1),
                  calc_MOI(POLYA), calc_MOI(PFPK2), calc_MOI(M2490))) %>%
  subset(MOI>0)
rownames(Angola_2021_data) <- Angola_2021_data$Sample_ID

# Nei's gene identity metric: sums of squares of allele frequencies
for (m in MARKERS) {
  Angola_2021_data[, paste0(m, "_diversity")] <- 
    1-apply(Angola_2021_data, 1, function(x) 
      expected_prop_IBS(as.numeric(x[["MOI"]]), calc_MOI(x[[m]])))
}

# ============================ PARSE CLINICAL METADATA ============================

Angola_2021_clinical <- read_xlsx("Dimbu_Angola_TES_metadata.xlsx") %>% as.data.frame
rownames(Angola_2021_clinical) <- Angola_2021_clinical$Codigo

Angola_2021_followup <- Angola_2021_clinical %>%
  dplyr::select(c("Codigo", paste0("Parasitemia d", c(0, 2,3,7,14,21,28,35,42)))) %>% 
  reshape2::melt(id="Codigo") %>% 
  mutate(variable=as.numeric(gsub("Parasitemia d", "", variable)), 
         value=as.numeric(value), 
         follow_up=ifelse(!is.na(value), variable, 0)) %>% 
  group_by(Codigo) %>% summarise(follow_up=max(follow_up)) %>% as.data.frame %>%
  transmute(pair=Codigo, Followup=follow_up,
            Province=Angola_2021_clinical[pair, "Provincia"],
            Drug=Angola_2021_clinical[pair, "Farmaco"],
            Dimbu_posterior=as.numeric(Angola_2021_clinical[pair, "Prob_Recr"]), 
            D0_parasitemia=Angola_2021_clinical[pair, "Parasitemia d0"],
            D0_temperature=Angola_2021_clinical[pair, "Temperatura d0"])

Angola_2021_followup$DX_parasitemia <- 
  as.numeric(apply(Angola_2021_followup, 1, function(x) {
    y <- Angola_2021_clinical[x[1], paste0("Parasitemia d", as.numeric(x[2]))]; 
    return(ifelse(is.null(y), NA, y))}))

Angola_2021_followup$DX_temperature <- 
  as.numeric(apply(Angola_2021_followup, 1, function(x) {
    y <- Angola_2021_clinical[x[1], paste0("Temperatura d", as.numeric(x[2]))]; 
    return(ifelse(is.null(y), NA, y))}))

Angola_2021_followup$DX_temperature <- 
  as.numeric(apply(Angola_2021_followup, 1, function(x) {
    y <- Angola_2021_clinical[x[1], paste0("Temperatura d", x[2])]; 
    return(ifelse(is.null(y), NA, y))}))

Angola_2021_followup$D0_MOI <- 
  as.numeric(apply(Angola_2021_followup, 1, function(x) {
    y <- Angola_2021_data %>% subset(Code==x[1] & Timepoint=="D0") 
    return(ifelse(is.null(y), NA, y[,"MOI"]))}))

Angola_2021_followup$DX_MOI <- 
  as.numeric(apply(Angola_2021_followup, 1, function(x) {
    y <- Angola_2021_data %>% subset(Code==x[1] & Timepoint!="D0") 
    return(ifelse(is.null(y), NA, y[,"MOI"]))}))

rownames(Angola_2021_followup) <- Angola_2021_followup$pair
Angola_2021_followup["ZQ21-103", "DX_parasitemia"] <- 
  as.numeric(Angola_2021_clinical["ZQ21-103", "Parasitemia dVNP1"])

Angola_2021_clin_metadata <- 
  merge(Angola_2021_followup[, c("pair", "D0_parasitemia", "DX_parasitemia", "D0_temperature", "DX_temperature")],
        Angola_2021_data[, c("Code", "Sample_ID", "Timepoint", "MOI")], by.x="pair", by.y="Code") %>% 
  transmute(Sample_ID=Sample_ID, Parasitemia=ifelse(Timepoint=="D0", D0_parasitemia, DX_parasitemia),
            Temperature=ifelse(Timepoint=="D0", D0_temperature, DX_temperature))
Angola_2021_data <- merge(Angola_2021_data, Angola_2021_clin_metadata)
rownames(Angola_2021_data) <- Angola_2021_data$Sample_ID

Angola_2021_pairs <- Angola_2021_data %>% 
  mutate(Timepoint=ifelse(grepl("D0", Timepoint), "D_0", "D_REC")) %>%
  reshape2::dcast(formula=Code~Timepoint, value.var="Sample_ID") %>%
  mutate(Dimbu_posterior=Angola_2021_data[D_0, "Plucinski_posterior"])

Angola_2021_pairs <- Angola_2021_pairs %>% subset(!is.na(D_REC))

# ============================ PARSE POSTERIORS ============================

Plucinski_classifications <- read_rds("../Plucinski_et_al_model/Full_dataset/Plucinski_classifications.rds")

Plucinski_posteriors <- lapply(TYPES, function(x) 
  Plucinski_classifications[[x]][["posterior_recrudesence"]] %>% 
    transmute(Code=isolate, !!x:=prob)) %>% Reduce(merge, .)

Angola_2021_pairs <- merge(Angola_2021_pairs, Plucinski_posteriors)
rownames(Angola_2021_pairs) <- Angola_2021_pairs$Code

write_rds(Angola_2021_data, "Dimbu_Angola_TES_data.rds")
write_rds(Angola_2021_pairs, "Dimbu_Angola_TES_pairs.rds")
write_rds(Angola_2021_followup, "Dimbu_Angola_TES_followup.rds")



# ============================ BASELINE METRICS ==============================

diversity_plot <- Angola_2021_data %>% 
  mutate(Type=ifelse(Timepoint=="D0", "Baseline", "Recurrent")) %>% 
  select(c("Parasitemia", "Temperature", "Type", paste0(MARKERS, "_diversity"))) %>% 
  reshape2::melt(id=c("Parasitemia", "Temperature", "Type")) %>% 
  mutate(variable=gsub("_diversity", "", variable)) %>%
  ggplot(aes(x=Parasitemia, y=value)) + 
  geom_smooth(method="lm", lwd=0.4, lty=2, alpha=0.25, col="#748ba1", fill="#97b0c7") +
  geom_point(aes(pch=Type, col=Type), size=2) + #, col="#545454") + 
  facet_wrap(vars(variable), nrow=2) + 
  scale_x_continuous(trans="log10") + coord_cartesian(ylim=c(0, 1)) +
  scale_shape_manual(values=c(1, 4)) + 
  scale_color_manual(values=c("#008080", "#6A5ACD")) +
  xlab("Parasitemia (parasites per μL)") +
  ylab("Within-host diversity*") +
  labs(caption=paste0("*Expected value of 1-Nei's gene identity metric, assuming M equifrequent clones ",
                      "(where M is the maximum cardinality across all loci)\nexhibit A distinct alleles ",
                      "(where A is the locus cardinality)")) + 
  theme_bw() + theme(strip.text = element_text(size=11),
                     legend.position = "none",
                     plot.caption = element_text(size=10.25, hjust=0.5))

MOI_covariates <- merge(Angola_2021_data, Angola_2021_clin_metadata) %>% 
  mutate(Fever=ifelse(Temperature>=37.5, "Temperature>=37.5", "Temperature<37.5"),
         Type=ifelse(Timepoint=="D0", "Baseline", "Recurrent")) %>% 
  subset(Parasitemia>0)

JT_overall <- 
  PMCMRplus::jonckheereTest(x=MOI_covariates$MOI, g=log10(MOI_covariates$Parasitemia), 
                            alternative = 'less')

MOI_covariates_plot <- ggplot() + 
  geom_boxplot(data=MOI_covariates, aes(x=Parasitemia, y=MOI, group=MOI), 
               outlier.size=0, size=0.25, alpha=0.25, fill="#97b0c7", col="#748ba1") + 
  geom_jitter(data=MOI_covariates, aes(x=Parasitemia, y=MOI, pch=Type, col=Type), 
              width=0, height=0.25, size=1.8, stroke=0.6) + 
  annotate("label", x=10^mean(log10(range(MOI_covariates$Parasitemia))), y=5.8,
           label=paste0("Jonckheere-Terpstra test for decreasing trend: p=", 
                        format(JT_overall[["p.value"]], digits=1, type="e"))) +
  scale_x_continuous(trans="log10") + scale_y_continuous(breaks=1:5) +
  scale_color_manual(values=c("#008080", "#6A5ACD")) +
  xlab("Parasitemia (parasites per μL)") +
  ylab("Multiplicity of infection\n(maximum cardinality\nacross markers)") +
  scale_shape_manual(values=c(1, 4)) + 
  theme_bw() + theme(strip.text = element_text(size=11), 
                     legend.position = "bottom",
                     legend.text = element_text(size=11))

pdf("Dimbu_Angola_MOI_covariates.pdf", height=8, width=10)
show(plot_grid(plot_grid(NULL, MOI_covariates_plot, NULL, nrow=1, rel_widths = c(1, 4, 1)), 
               diversity_plot, ncol=1, rel_heights = c(2, 3), 
               align="h", axis="lr", labels=c("(A)", "(B)")))
dev.off()



