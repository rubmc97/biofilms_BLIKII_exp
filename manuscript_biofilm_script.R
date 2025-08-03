library(phyloseq)
library(decontam)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ANCOMBC)
library(tidyr)
library(NetCoMi)
library(stringr)
library(dada2)
library(plotly)
library(knitr)
library(Biostrings)

#Setting the path where the files are after the fastp processing to remove the adapters
setwd("/home/comi/ruben.martinez/biofilm_samples/")
path = "/home/comi/ruben.martinez/biofilm_samples/"
list.files(path)

#Detection of forward and reverse and setting the sample names
fnFs = sort(list.files(path, pattern="_R1_fastp.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2_fastp.fastq.gz", full.names = TRUE))

sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

forwplot = ggplotly(plotQualityProfile(fnFs[1:13], n = 1e+05, aggregate = TRUE) +
                      geom_hline(yintercept=c(20,25,30), color=c("red","blue","green"), size=0.5),
                    width =750)
forwplot

revqplot = ggplotly(plotQualityProfile(fnRs[1:13], n = 1e+05, aggregate = T) + 
                      geom_hline(yintercept=c(20,25,30),
                                 color=c("red","blue","green"),
                                 size=0.2)) %>% ggplotly(800)
revqplot

# Place filtered files in filtered/ subdirectory
filtFs = file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(263,201),
                    maxN=0, maxEE=c(3,4), truncQ=2, rm.phix=TRUE, trimLeft=c(19,20),
                    compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

errF = learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

errR = learnErrors(filtRs, multithread=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Dereplication step
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)
names(derepFs) = sample.names
names(derepRs) = sample.names

#Sample Inference
dadaFs = dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]

dadaRs = dada(filtRs, err=errR, multithread=TRUE)
dadaRs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab2 = seqtab[,nchar(colnames(seqtab)) %in% 285:335] #determine depending on the previous step

#Remove chimeric sequences
seqtab.nochim = removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, file = "/home/comi/ruben.martinez/seqtab.nochim.rds")
save(seqtab.nochim, file = "/home/comi/ruben.martinez/seqtab.nochim.RData")

getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample.names
head(track)

write.csv(track, "/home/comi/ruben.martinez/track_table_biofilms_0624.csv")
#Assign taxonomy using DECIPHER and the trained classifier from the SILVA database

library(DECIPHER); packageVersion("DECIPHER")
dna = DNAStringSet(getSequences(seqtab.nochim))
load("SILVA_SSU_r132_March2018.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid = t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) = ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa = taxid

#Handoff to phyloseq

sd_biofilms = read.csv("metada_biofilms.csv")

phyloseq_biofilms = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                             sample_data(sd_biofilms), 
                             tax_table(taxa))

dna = Biostrings::DNAStringSet(taxa_names(phyloseq_biofilms))
names(dna) = taxa_names(phyloseq_biofilms)
ps = merge_phyloseq(phyloseq_biofilms, dna)
taxa_names(phyloseq_biofilms) = paste0("ASV", seq(ntaxa(phyloseq_biofilms)))
phyloseq_biofilms


  load("phyloseq_biofilms_no_filt_150425.RData")
  
track = read.csv("track_table_biofilms_0425.csv", row.names = 1)
track$SampleID = row.names(track)
sd = as.matrix.data.frame(sample_data(phyloseq_biofilms))

track.sd = merge(track, sd, by.y = "SampleID")

write.csv(track.sd, "TableS1_supp_track.csv")

sample_data(phyloseq_biofilms)$is.neg = sample_data(phyloseq_biofilms)$site == "negative"

contamdf.prev05= isContaminant(phyloseq_biofilms, method="prevalence", neg="is.neg", threshold=0.05)

table(contamdf.prev05$contaminant)

phyloseq.biofilms.filt = prune_taxa(!contamdf.prev05$contaminant, phyloseq_biofilms) #remove contaminants from the negatives

#physeq_blik_2_filt = phy_blik_2_filt %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Class) ) #remove sequences labelled as mitochondria or chloroplast

phyloseq.biofilms.filt = phyloseq.biofilms.filt %>%
  subset_taxa(
    !grepl("Chloroplast", family, ignore.case = TRUE) &
      !grepl("mitochondria", family, ignore.case = TRUE) &
      !grepl("Chloroplast", order, ignore.case = TRUE) &
      !grepl("mitochondria", order, ignore.case = TRUE) &
      !grepl("Chloroplast", class, ignore.case = TRUE) &
      !grepl("mitochondria", class, ignore.case = TRUE) &
      !grepl("Chloroplast", phylum, ignore.case = TRUE) &
      !grepl("mitochondria", phylum, ignore.case = TRUE)
  )

phyloseq.biofilms.filt = subset_samples(phyloseq.biofilms.filt, is.neg !=TRUE)
sample_data(phyloseq.biofilms.filt)$is.neg = NULL

# First, extract the tax_table
tax = tax_table(phyloseq.biofilms.filt)

# Find ASVs where domain is not NA
keep_taxa = as.vector(!is.na(tax[, "domain"]))

# Subset phyloseq to keep only those ASVs
phyloseq.biofilms.filt = prune_taxa(keep_taxa, phyloseq.biofilms.filt)
phyloseq.biofilms.filt = prune_taxa(taxa_sums(phyloseq.biofilms.filt) > 0, phyloseq.biofilms.filt)

save(phyloseq.biofilms.filt,
     file = paste("phyloseq_biofilms_filt_16042025",
                  ".RData",
                  sep = ""))

otu_table(phyloseq.biofilms.filt) = otu_table(
  t(otu_table(phyloseq.biofilms.filt)),
  taxa_are_rows = TRUE
)

#Plot rarefaction curves
tab = otu_table(phyloseq.biofilms.filt)
class(tab) = "matrix"
tab = t(tab)
rare = rarecurve(tab, step=10000, lwd=2, ylab="ASVs",  label=F)

phy.blik2.css = microbiomeMarker::normalize(phyloseq.biofilms.filt, method = "CSS")
phy.blik2.css = prune_taxa(taxa_sums(phy.blik2.css) > 0, phy.blik2.css)

save(phy.blik2.css,
     file = paste("phy_biofilms_norm_16042025",
                  ".RData",
                  sep = ""))

load("phy_biofilms_norm_16042025.RData")
##calculate the total ASVs number per sample and plot the results
asvs.number = microbiome::alpha(phy.blik2.css, "observed")
asvs.number$SampleID = row.names(asvs.number)
sd.blik2 = sample_data(phy.blik2.css)
sd.blik2 = as.matrix.data.frame(sd.blik2)

df.asvs.number = merge(asvs.number, sd.blik2, by.y = "SampleID")

alpha.ott = subset.data.frame(df.asvs.number, stream == "otterbach")
# Load the dplyr package
library(dplyr)

# Assuming your data frame is named df
df_avg = alpha.ott %>%
  group_by(site, st, sample_type) %>%
  summarise(avg_observed = mean(observed, na.rm = TRUE))

as.alpha.ott = subset.data.frame(alpha.ott, sample_type =="AS")
ns.alpha.ott = subset.data.frame(alpha.ott, sample_type =="NS")

#plot for the Otterbach dataset
custom_order = c("EO", "IO", "FO")

colors = c("#999999", "#009E73", "#E69F00")

as.plot.alpha.ott = ggplot(as.alpha.ott, aes(x = factor(site, levels = custom_order), y = observed, fill = land_use)) +
  geom_boxplot(alpha = 0.5, width = 0.5) +
  geom_jitter(size = 1, shape=16, alpha = 0.6) +
  facet_grid( ~ st) +
  scale_fill_manual(values = colors) +  # Set custom colors
  ggpubr::theme_pubclean() + labs(x = "Site", y = "Observed ASVs") +
  ggpubr::stat_compare_means(comparisons = list(c("EO", "IO"), c("IO", "FO"), c("EO", "FO")), label = "p.signif", hide.ns = TRUE, method = "wilcox.test", p.adjust.methods = "fdr") +
  theme(
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    axis.text = element_text(color = "black")   # Ensure axis text is visible
  ) +
  scale_y_continuous(limits = c(0, 2000), breaks = c(500, 1000, 1500, 2000))
as.plot.alpha.ott

ns.plot.alpha.ott = ggplot(ns.alpha.ott, aes(x = factor(site, levels = custom_order), y = observed, fill = land_use)) +
  geom_boxplot(alpha = 0.5, width = 0.5) +
  geom_jitter(size = 1, shape=16, alpha = 0.6) +
  facet_grid( ~ st) +
  scale_fill_manual(values = colors) +  # Set custom colors
  ggpubr::theme_pubclean() + labs(x = "Site", y = "Observed ASVs") +
  ggpubr::stat_compare_means(comparisons = list(c("EO", "IO"), c("IO", "FO"), c("EO", "FO")), label = "p.signif", hide.ns = TRUE, method = "wilcox.test", p.adjust.methods = "fdr") +
  theme(
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    axis.text = element_text(color = "black")   # Ensure axis text is visible
  ) +
  scale_y_continuous(limits = c(0, 2000))
ns.plot.alpha.ott

as.ns.ott = ggpubr::ggarrange(as.plot.alpha.ott, ns.plot.alpha.ott, ncol = 2, labels=c("A"), common.legend = T)
as.ns.ott

#observed ASVs plot of the Perlenbach dataset
alpha.pb = subset.data.frame(df.asvs.number, stream == "perlenbach")
as.alpha.pb = subset.data.frame(alpha.pb, sample_type == "AS")
ns.alpha.pb = subset.data.frame(alpha.pb, sample_type == "NS")
#plot for the Perlenbach dataset
colors_pb = c("#999999", "#009E73")
facet_order = c("AS", "NS")
custom_order_pb = c("EP", "FO")

alpha.pb$sample_type = factor(alpha.pb$sample_type, levels = facet_order)

as.plot.alpha.pb = ggplot(as.alpha.pb, aes(x = factor(site, levels = custom_order_pb), y = observed, fill = land_use)) +
  geom_boxplot(alpha = 0.5, width = 0.35) +
  geom_jitter(size = 1, shape=16, alpha = 0.6) +
  facet_grid(~ st) +
  labs(x = "Site", y = "Observed ASVs") +
  scale_fill_manual(values = colors_pb) +  # Set custom colors
  ggpubr::theme_pubclean() +
  ggpubr::stat_compare_means(comparisons = list(c("EP", "FO")), label = "p.signif", hide.ns = TRUE, method = "wilcox.test", p.adjust.methods = "fdr") +
  theme(
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    axis.text = element_text(color = "black")   # Ensure axis text is visible
  ) +
  scale_y_continuous(limits = c(250, 1250), breaks = c(250, 500, 750, 1000, 1250))
as.plot.alpha.pb

ns.plot.alpha.pb = ggplot(ns.alpha.pb, aes(x = factor(site, levels = custom_order_pb), y = observed, fill = land_use)) +
  geom_boxplot(alpha = 0.5, width = 0.35) +
  geom_jitter(size = 1, shape=16, alpha = 0.6) +
  facet_grid(~ st) +
  labs(x = "Site", y = "Observed ASVs") +
  scale_fill_manual(values = colors_pb) +  # Set custom colors
  ggpubr::theme_pubclean() +
  ggpubr::stat_compare_means(comparisons = list(c("EP", "FO")), label = "p.signif", hide.ns = TRUE, method = "wilcox.test", p.adjust.methods = "fdr") +
  theme(
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    axis.text = element_text(color = "black")   # Ensure axis text is visible
  ) +
  scale_y_continuous(limits = c(250, 1250), breaks = c(250, 500, 750, 1000, 1250))
ns.plot.alpha.pb

as.ns.pb = ggpubr::ggarrange(as.plot.alpha.pb, ns.plot.alpha.pb, ncol = 2, common.legend = T)
as.ns.ott.pb = ggpubr::ggarrange(as.ns.ott, as.ns.pb, nrow = 2, common.legend = T)
as.ns.ott.pb

#Ordination with PCoA and dissimilarity with Bray-Curtis
blik.ott = subset_samples(phy.blik2.css, stream == "otterbach")
as.ott = subset_samples(blik.ott, sample_type == "AS")
ns.ott = subset_samples(blik.ott, sample_type == "NS")
as.ott = prune_taxa(taxa_sums(as.ott) > 0, as.ott)
ns.ott = prune_taxa(taxa_sums(ns.ott) > 0, ns.ott)

bray.as.ott = phyloseq::distance(as.ott, method = "bray")
pcoa.res = ape::pcoa(bray.as.ott)
pcoa.as.ott.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                            Axis1 = pcoa.res$vectors[, 1],
                            Axis2 = pcoa.res$vectors[, 2])

sd.ott.as = sample_data(as.ott) %>% as.matrix.data.frame()
write.csv(sd.ott.as, "sd_as_ott.csv")
sd.ott.as = read.csv("sd_as_ott.csv", row.names = 1)


pcoa.df = merge(pcoa.as.ott.df, sd.ott.as, by.x = "SampleID", by.y = "row.names")
pcoa.df$SampleID.y = NULL
pcoa.df$sample_type = NULL
pcoa.df$SampleID = NULL

nut.ott.as = read.csv("nuts_as_ott.csv", row.names = 1)
all(rownames(sample_data(as.ott)) %in% rownames(nut.ott.as))

env.fit = envfit(pcoa.res$vectors, nut.ott.as, permutations = 999)
# Extract p-values
pvals = env.fit$vectors$pvals

# Get names of variables with p < 0.05
padj = p.adjust(pvals, method = "fdr")
significant.vars = names(padj[padj < 0.05])

# Filter the vectors
significant.vectors = env.fit$vectors$arrows[significant.vars, ]

vec_lengths = env.fit$vectors$r 
env.vectors = as.data.frame(env.fit$vectors$arrows)
env.vectors$Variable = rownames(env.vectors)
env.vectors$Length = env.fit$vectors$r
env.vectors$x_end = env.vectors$Axis.1 * env.vectors$Length * 0.5
env.vectors$y_end = env.vectors$Axis.2 * env.vectors$Length * 0.5


lu.colors = c("intensive" = "#E69F00", "forest" = "#009E73", "extensive" = "#999999")
land.use = sample_data(as.ott)$land_use
st = sample_data(as.ott)$st

pcoa.as.ott = ggplot(pcoa.df, aes(x = Axis1, y = Axis2)) +  # Adjust "LandUse" as needed
  geom_point(aes(color = factor(land_use), shape = factor(st)), size = 3, alpha = 0.6) + scale_color_manual(values = lu.colors) +
  geom_segment(data = env.vectors,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = env.vectors,
                  aes(x = x_end, y = y_end, label = Variable),
                  size = 4) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA AS") + scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) + theme_bw()
pcoa.as.ott

bray.ns.ott = phyloseq::distance(ns.ott, method = "bray")

pcoa.res = ape::pcoa(bray.ns.ott)

#write.csv(as.matrix.data.frame(sample_data(ns.ott)), "sd_ns_ott.csv")
pcoa.ns.ott.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                            Axis1 = pcoa.res$vectors[, 1],
                            Axis2 = pcoa.res$vectors[, 2])

sd.ott.ns = read.csv("sd_ns_ott.csv", row.names = 1)

pcoa.df = merge(pcoa.ns.ott.df, sd.ott.ns, by.x = "SampleID", by.y = "row.names")
pcoa.df$river = NULL
pcoa.df$sample_type = NULL

row.names(pcoa.df) = pcoa.df$SampleID

nut.ott.ns = read.csv("nuts_ns_ott.csv", row.names = 1)
all(rownames(sample_data(ns.ott)) %in% rownames(nut.ott.ns))

env.fit = envfit(pcoa.res$vectors, nut.ott.ns, permutations = 999)
pvals = env.fit$vectors$pvals
# Get names of variables with p < 0.05
padj = p.adjust(pvals, method = "fdr")
significant.vars = names(padj[padj < 0.05])

# Filter the vectors
significant.vectors = env.fit$vectors$arrows[significant.vars, ]

vec_lengths = env.fit$vectors$r 
env.vectors = as.data.frame(env.fit$vectors$arrows)
env.vectors$Variable = rownames(env.vectors)
env.vectors$Length = env.fit$vectors$r
env.vectors$x_end = env.vectors$Axis.1 * env.vectors$Length * 0.5
env.vectors$y_end = env.vectors$Axis.2 * env.vectors$Length * 0.5


lu.colors = c("intensive" = "#E69F00", "forest" = "#009E73", "extensive" = "#999999")
land.use = sample_data(ns.ott)$land_use
st = sample_data(ns.ott)$st

pcoa.ns.ott = ggplot(pcoa.df, aes(x = Axis1, y = Axis2)) +  # Adjust "LandUse" as needed
  geom_point(aes(color = factor(land_use), shape = factor(st)), size = 3, alpha = 0.6) + scale_color_manual(values = lu.colors) +
  geom_segment(data = env.vectors,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = env.vectors,
                  aes(x = x_end, y = y_end, label = Variable),
                  size = 4) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA NS") + scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) + theme_bw()
pcoa.ns.ott

blik.pb = subset_samples(phy.blik2.css, stream == "perlenbach")
as.pb = subset_samples(blik.pb, sample_type == "AS")
ns.pb = subset_samples(blik.pb, sample_type == "NS")
as.pb = prune_taxa(taxa_sums(as.pb) > 0, as.pb)
ns.pb = prune_taxa(taxa_sums(ns.pb) > 0, ns.pb)

bray.as.pb = phyloseq::distance(as.pb, method = "bray")
pcoa.res = ape::pcoa(bray.as.pb)

pcoa.as.pb.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                           Axis1 = pcoa.res$vectors[, 1],
                           Axis2 = pcoa.res$vectors[, 2])
#write.csv(as.matrix.data.frame(sample_data(as.pb)), "sd_pb_as.csv")
sd.pb.as = read.csv("sd_pb_as.csv", row.names = 1)

pcoa.df = merge(pcoa.as.pb.df, sd.pb.as, by.x = "SampleID", by.y = "row.names")
row.names(pcoa.df) = pcoa.df$SampleID

nut.pb.as = read.csv("nuts_as_pb.csv", row.names = 1)
all(rownames(sample_data(as.pb)) %in% rownames(nut.pb.as))

env.fit = envfit(pcoa.res$vectors, nut.pb.as, permutations = 999)
# Extract p-values
pvals = env.fit$vectors$pvals

# Get names of variables with p < 0.05
padj = p.adjust(pvals, method = "fdr")
significant.vars = names(padj[padj < 0.05])

# Filter the vectors
significant.vectors = env.fit$vectors$arrows[significant.vars, ]

vec_lengths = env.fit$vectors$r 
env.vectors = as.data.frame(env.fit$vectors$arrows)
env.vectors$Variable = rownames(env.vectors)
env.vectors$Length = env.fit$vectors$r
env.vectors$x_end = env.vectors$Axis.1 * env.vectors$Length * 0.5
env.vectors$y_end = env.vectors$Axis.2 * env.vectors$Length * 0.5


lu.colors = c("forest" = "#009E73", "extensive" = "#999999")
land.use = sample_data(as.pb)$land_use
st = sample_data(as.pb)$st

pcoa.as.pb = ggplot(pcoa.df, aes(x = Axis1, y = Axis2)) +  # Adjust "LandUse" as needed
  geom_point(aes(color = factor(land_use), shape = factor(st)), size = 3, alpha = 0.6) + scale_color_manual(values = lu.colors) +
  geom_segment(data = env.vectors,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = env.vectors,
                  aes(x = x_end, y = y_end, label = Variable),
                  size = 4) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA AS") + scale_x_continuous(limits = c(-0.55, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) + theme_bw()
pcoa.as.pb

bray.ns.pb = phyloseq::distance(ns.pb, method = "bray")
pcoa.res = ape::pcoa(bray.ns.pb)

pcoa.ns.pb.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                           Axis1 = pcoa.res$vectors[, 1],
                           Axis2 = pcoa.res$vectors[, 2])
#write.csv(as.matrix.data.frame(sample_data(ns.pb)), "sd_pb_ns.csv")
sd.pb.ns = read.csv("sd_pb_ns.csv", row.names = 1)

pcoa.df = merge(pcoa.ns.pb.df, sd.pb.ns, by.x = "SampleID", by.y = "row.names")
row.names(pcoa.df) = pcoa.df$SampleID

nut.pb.ns = read.csv("nuts_ns_pb.csv", row.names = 1)
all(rownames(sample_data(ns.pb)) %in% rownames(nut.pb.ns))

env.fit = envfit(pcoa.res$vectors, nut.pb.ns, permutations = 999)
# Extract p-values
pvals = env.fit$vectors$pvals

# Get names of variables with p < 0.05
padj = p.adjust(pvals, method = "fdr")
significant.vars = names(padj[padj < 0.05])

# Filter the vectors
significant.vectors = env.fit$vectors$arrows[significant.vars, ]

vec_lengths = env.fit$vectors$r 
env.vectors = as.data.frame(env.fit$vectors$arrows)
env.vectors$Variable = rownames(env.vectors)
env.vectors$Length = env.fit$vectors$r
env.vectors$x_end = env.vectors$Axis.1 * env.vectors$Length * 0.5
env.vectors$y_end = env.vectors$Axis.2 * env.vectors$Length * 0.5

lu.colors = c("forest" = "#009E73", "extensive" = "#999999")

pcoa.ns.pb = ggplot(pcoa.df, aes(x = Axis1, y = Axis2)) +  # Adjust "LandUse" as needed
  geom_point(aes(color = factor(land_use), shape = factor(st)), size = 3, alpha = 0.6) + scale_color_manual(values = lu.colors) +
  geom_segment(data = env.vectors,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = env.vectors,
                  aes(x = x_end, y = y_end, label = Variable),
                  size = 4) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA NS") + theme_bw() + scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25))
pcoa.ns.pb

all.pcoa = ggpubr::ggarrange(pcoa.as.ott, pcoa.ns.ott, pcoa.as.pb, pcoa.ns.pb, ncol = 2, nrow = 2, common.legend = T, labels = "B")
all.pcoa

#PERMANOVA FOR AS OTT
adonis2(bray.as.ott ~land_use, data = sd.ott.as)
adonis2(bray.as.ott ~st, data = sd.ott.as)
adonis2(bray.as.ott ~land_use * st, data = sd.ott.as)

#PERMANOVA FOR NS OTT
adonis2(bray.ns.ott ~land_use, data = sd.ott.ns)
adonis2(bray.ns.ott ~st, data = sd.ott.ns)
adonis2(bray.ns.ott ~land_use * st, data = sd.ott.ns)

#PERMANOVA FOR AS PB
adonis2(bray.as.pb ~land_use, data = sd.pb.as)
adonis2(bray.as.pb ~st, data = sd.pb.as)
adonis2(bray.as.pb ~land_use * st, data = sd.pb.as)

#PERMANOVA FOR NS PB
adonis2(bray.ns.pb ~land_use, data = sd.pb.ns)
adonis2(bray.ns.pb ~st, data = sd.pb.ns)
adonis2(bray.ns.pb ~land_use * st, data = sd.pb.ns)

#Subsetting samples for the ANCOMBC2 anaylses
as.ott.t1 = subset_samples(as.ott, st == "T1")
as.ott.t2 = subset_samples(as.ott, st == "T2")

ns.ott.t1 = subset_samples(ns.ott, st == "T1")
ns.ott.t2 = subset_samples(ns.ott, st == "T2")

as.pb.t1 = subset_samples(as.pb, st == "T1")
as.pb.t1 = prune_taxa(taxa_sums(as.pb.t1) > 0, as.pb.t1)
as.pb.t2 = subset_samples(as.pb, st == "T2")
as.pb.t2 = prune_taxa(taxa_sums(as.pb.t2) > 0, as.pb.t2)

ns.pb.t1 = subset_samples(ns.pb, st == "T1")
ns.pb.t1 = prune_taxa(taxa_sums(ns.pb.t1) > 0, ns.pb.t1)
ns.pb.t2 = subset_samples(ns.pb, st == "T2")
ns.pb.t2 = prune_taxa(taxa_sums(ns.pb.t2) > 0, ns.pb.t2)

#Running ANCOMBC2
load("phyloseq_biofilms_filt_16042025.RData") #load the phyloseq object prior to CSS normalization
blik.ott = subset_samples(phyloseq.biofilms.filt, stream == "otterbach")
as.ott = subset_samples(blik.ott, sample_type == "AS")
ns.ott = subset_samples(blik.ott, sample_type == "NS")
as.ott = prune_taxa(taxa_sums(as.ott) > 0, as.ott)
ns.ott = prune_taxa(taxa_sums(ns.ott) > 0, ns.ott)

sd.as.ott = phyloseq::sample_data(as.ott)

sample_data(as.ott)$land_use_st = interaction(sample_data(as.ott)$land_use, sample_data(as.ott)$st)
# Make sure the combined factor is a factor
sample_data(as.ott)$land_use_st = factor(sample_data(as.ott)$land_use_st)

# Set 'extensive_T2' as the reference level
sample_data(as.ott)$land_use_st = relevel(sample_data(as.ott)$land_use_st, ref = "forest.T1")
sample_data(as.ott)$land_use_st = relevel(sample_data(as.ott)$land_use_st, ref = "forest.T2")

ancombc2.as.t1 = ancombc2(as.ott, p_adj_method = "fdr", fix_formula = "land_use_st", group = "land_use_st", tax_level = NULL,
                          alpha = 0.05, verbose = F, dunnet = T, pairwise = F)

ancombc2.as.t2 = ancombc2(as.ott, p_adj_method = "fdr", fix_formula = "land_use_st", group = "land_use_st", tax_level = NULL,
                          alpha = 0.05, verbose = F, dunnet = T, pairwise = F)

res.ancombc2.as.t1 = ancombc2.as.t1$res
res.ancombc2.as.t2 = ancombc2.as.t2$res

res.ancombc2.as.t1$ASV = res.ancombc2.as.t1$taxon
res.ancombc2.as.t2$ASV = res.ancombc2.as.t2$taxon

res.ancombc2.as.t1 = res.ancombc2.as.t1[, c(2, 4, 11, 13, 23), drop = FALSE] #select columns of interest

res.ancombc2.as.t1_1 = res.ancombc2.as.t1 %>%
  filter((p_land_use_stextensive.T1 < 0.05 &
            abs(lfc_land_use_stextensive.T1) > 0.5) |
           (p_land_use_stintensive.T1 < 0.05 &
              abs(lfc_land_use_stintensive.T1) > 0.5))

res.ancombc2.as.t2_1 = res.ancombc2.as.t2 %>%
  filter((p_land_use_stextensive.T2 < 0.05 &
            abs(lfc_land_use_stextensive.T2) > 0.5) |
           (p_land_use_stintensive.T2 < 0.05 &
              abs(lfc_land_use_stintensive.T2) > 0.5))

res.ancombc2.as.t1.mut = res.ancombc2.as.t1 %>%
  pivot_longer(cols = c(lfc_land_useforest, lfc_land_useintensive_land_useforest), 
               names_to = "lfc_source", values_to = "lfc_value") %>%
  mutate(lfc_value = ifelse(lfc_source == "lfc_land_useforest", -lfc_value, lfc_value))

res.ancombc2.as.t1.mut.long = res.ancombc2.as.t1.mut %>%
  pivot_longer(
    cols = starts_with("p_land_use"),
    names_to = "contrast",
    values_to = "p_value"
  )

#filtering the results for p value and lfc
res.ancombc2.as.t1.mut.long = res.ancombc2.as.t1.mut.long %>%
  filter(abs(p_value) < 0.05)

res.ancombc2.as.t1.mut.long = res.ancombc2.as.t1.mut.long %>%
  filter(abs(lfc_value) > 0.5)

res.ancombc2.as.t1.mut.long.filt = res.ancombc2.as.t1.mut.long %>%
  filter(
    str_remove(contrast, "^[^_]+_") == str_remove(lfc_source, "^[^_]+_")
  )

tt.as.ott = as.data.frame(tax_table(as.ott))
tt.as.ott$ASV = row.names(tt.as.ott)

df.res.ancombc2.as.t1 = merge(res.ancombc2.as.t1_1, tt.as.ott, by.y = "ASV")

ancombc2.as.t2 = ancombc2(as.ott.t2, p_adj_method = "fdr", fix_formula = "land_use", group = "land_use", tax_level = NULL, alpha = 0.05, verbose = F, dunnet = F, pairwise = T)

res.ancombc2.as.t2 = ancombc2.as.t2$res_pair

res.ancombc2.as.t2$ASV = res.ancombc2.as.t2$taxon

res.ancombc2.as.t2 = res.ancombc2.as.t2[, c(2, 4, 11, 13, 14, 16, 23), drop = FALSE] #select columns of interest

res.ancombc2.as.t2.mut = res.ancombc2.as.t2 %>%
  pivot_longer(cols = c(lfc_land_useforest, lfc_land_useintensive_land_useforest), 
               names_to = "lfc_source", values_to = "lfc_value") %>%
  mutate(lfc_value = ifelse(lfc_source == "lfc_land_useforest", -lfc_value, lfc_value))

res.ancombc2.as.t2.mut.long = res.ancombc2.as.t2.mut %>%
  pivot_longer(
    cols = starts_with("q_land_use"),
    names_to = "contrast",
    values_to = "q_value"
  )

#filtering the results for p value and lfc
res.ancombc2.as.t2.mut.long = res.ancombc2.as.t2.mut.long %>%
  filter(abs(p_value) < 0.05)

res.ancombc2.as.t2.mut.long = res.ancombc2.as.t2.mut.long %>%
  filter(abs(lfc_value) > 0.5)

res.ancombc2.as.t2.mut.long.filt = res.ancombc2.as.t2.mut.long %>%
  filter(
    str_remove(contrast, "^[^_]+_") == str_remove(lfc_source, "^[^_]+_")
  )

tt.as.ott.t2 = as.data.frame(tax_table(as.ott.t2))
tt.as.ott.t2$ASV = row.names(tt.as.ott.t2)

res.ancombc2.as.ott.t2.df = merge(res.ancombc2.as.t2.mut.long.filt, tt.as.ott.t2, by.y = "ASV")

ancombc2.ns.t1 = ancombc2(ns.ott.t1, p_adj_method = "fdr", fix_formula = "land_use", group = "land_use", tax_level = NULL, alpha = 0.05, verbose = F, dunnet = F, pairwise = T)

res.ancombc2.ns.t1 = ancombc2.ns.t1$res_pair

res.ancombc2.ns.t1$ns. = res.ancombc2.ns.t1$taxon

res.ancombc2.ns.t1 = res.ancombc2.ns.t1[, c(2, 4, 11, 13, 23), drop = FALSE] #select columns of interest

res.ancombc2.ns.t1.mut = res.ancombc2.ns.t1 %>%
  pivot_longer(cols = c(lfc_land_useforest, lfc_land_useintensive_land_useforest), 
               names_to = "lfc_source", values_to = "lfc_value") %>%
  mutate(lfc_value = ifelse(lfc_source == "lfc_land_useforest", -lfc_value, lfc_value))

res.ancombc2.ns.t1.mut.long = res.ancombc2.ns.t1.mut %>%
  pivot_longer(
    cols = starts_with("p_land_use"),
    names_to = "contrast",
    values_to = "p_value"
  )

#filtering the results for p value and lfc
res.ancombc2.ns.t1.mut.long = res.ancombc2.ns.t1.mut.long %>%
  filter(abs(p_value) < 0.05)

res.ancombc2.ns.t1.mut.long = res.ancombc2.ns.t1.mut.long %>%
  filter(abs(lfc_value) > 0.5)

res.ancombc2.ns.t1.mut.long.filt = res.ancombc2.ns.t1.mut.long %>%
  filter(
    str_remove(contrns., "^[^_]+_") == str_remove(lfc_source, "^[^_]+_")
  )

ancombc2.ns.t2 = ancombc2(ns.ott.t2, p_adj_method = "fdr", fix_formula = "land_use", group = "land_use", tax_level = NULL, alpha = 0.05, verbose = F, dunnet = F, pairwise = T)

res.ancombc2.ns.t2 = ancombc2.ns.t2$res_pair

res.ancombc2.ns.t2$ns. = res.ancombc2.ns.t2$taxon

res.ancombc2.ns.t2 = res.ancombc2.ns.t2[, c(2, 4, 11, 13, 23), drop = FALSE] #select columns of interest

res.ancombc2.ns.t2.mut = res.ancombc2.ns.t2 %>%
  pivot_longer(cols = c(lfc_land_useforest, lfc_land_useintensive_land_useforest), 
               names_to = "lfc_source", values_to = "lfc_value") %>%
  mutate(lfc_value = ifelse(lfc_source == "lfc_land_useforest", -lfc_value, lfc_value))

res.ancombc2.ns.t2.mut.long = res.ancombc2.ns.t2.mut %>%
  pivot_longer(
    cols = starts_with("p_land_use"),
    names_to = "contrast",
    values_to = "p_value"
  )

#filtering the results for p value and lfc
res.ancombc2.ns.t2.mut.long = res.ancombc2.ns.t2.mut.long %>%
  filter(abs(p_value) < 0.05)

res.ancombc2.ns.t2.mut.long = res.ancombc2.ns.t2.mut.long %>%
  filter(abs(lfc_value) > 0.5)

res.ancombc2.ns.t2.mut.long.filt = res.ancombc2.ns.t2.mut.long %>%
  filter(
    str_remove(contrns., "^[^_]+_") == str_remove(lfc_source, "^[^_]+_")
  )

write.csv(res.ancombc2.as.t1.mut.long.filt, "ancombc2_AS_T1.csv")
write.csv(res.ancombc2.as.t2.mut.long.filt, "ancombc2_AS_T2.csv")
write.csv(res.ancombc2.ns.t1.mut.long.filt, "ancombc2_NS_T1.csv")
write.csv(res.ancombc2.ns.t2.mut.long.filt, "ancombc2_NS_T2.csv")

#Perlenbach samples

ancombc2.as.pb <- ancombc2(
  as.pb,
  p_adj_method = "fdr",
  fix_formula = "land_use_st",
  tax_level = NULL,
  alpha = 0.05,
  verbose = TRUE,
  dunnet = FALSE,
  pairwise = FALSE,
  ref
)

ancombc2.ns.pb.t2 = ancombc2(ns.pb.t2, p_adj_method = "fdr", fix_formula = "land_use", tax_level = NULL,
                             alpha = 0.05, verbose = T, dunnet = F, pairwise = F)


res.ancombc2.as.pb.T1 = ancombc2.as.pb.T1$res

res.ancombc2.as.pb.T2 = ancombc2.as.pb.T2$res

res.ancombc2.as.pb.T1$ASV = res.ancombc2.as.pb.T1$taxon

res.ancombc2.as.pb.T2$ASV = res.ancombc2.as.pb.T2$taxon

res.ancombc2.as.pb = res.ancombc2.as.pb[, c(3, 15, 23), drop = FALSE] #select columns of interest

#filtering the results for q value, ss and lfc
res.ancombc2.as.pb.t1.filt = res.ancombc2.as.pb %>%
  filter(p_land_use_stforest.T1 < 0.05) %>%
  filter(passed_ss_land_use_stforest.T1 == TRUE) %>%
  filter(abs(lfc_land_use_stforest.T1) > 0.5)

res.ancombc2.as.pb.t2.filt = res.ancombc2.as.pb.T2 %>%
  filter(p_land_use_stforest.T2 < 0.05) %>%
  filter(passed_ss_land_use_stforest.T2 == TRUE) %>%
  filter(abs(lfc_land_use_stforest.T2) > 0.5)

res.ancombc2.as.pb.t1.filt = res.ancombc2.as.pb.t1.filt[, c(3, 15, 30), drop = FALSE] #select columns of interest

res.ancombc2.as.pb.t2.filt = res.ancombc2.as.pb.t2.filt[, c(5, 17, 30), drop = FALSE] #select columns of interest

tt.as.pb = as.data.frame(tax_table(as.pb))
tt.as.pb$ASV = row.names(tt.as.pb)

res.ancombc2.as.pb.t1.filt.df = merge(res.ancombc2.as.pb.t1.filt, tt.as.pb, by.y = "ASV")
res.ancombc2.as.pb.t2.filt.df = merge(res.ancombc2.as.pb.t2.filt, tt.as.pb, by.y = "ASV")

write.csv(res.ancombc2.as.pb.t1.filt.df, "res_ancombc2_perlenbach_T1.csv")
write.csv(res.ancombc2.as.pb.t2.filt.df, "res_ancombc2_perlenbach_T2.csv")

phyla.colors <- c(
  "Acidobacteriota"   = "#556B2F",  # dark olive green
  "Actinomycetota"  = "#8B7500",  # dark golden brown
  "Bacteroidota"      = "#305A91",  # desaturated dark blue
  "Deinococcota"      = "#8B2323",  # dark brick red
  "Gemmatimonadota"   = "#2E8B57",  # deep sea green
  "Myxococcota"       = "#5E3C99",  # muted dark purple
  "Nitrospirota"      = "#7B3F61",  # dark muted pink-purple
  "Pseudomonadota"    = "#8B4513"   # saddle brown
)

setwd("/Dokumente und Einstellungen/RubenMartinezCuesta/Desktop/periphyton_subset_brief_report/")
anbc2.as.t1 = read.csv("res_ancombc2_AS_T1_wt.csv")
anbc2.as.t2 = read.csv("res_ancombc2_AS_T2_wt.csv")
anbc2.ns.t1 = read.csv("res_ancombc2_NS_T1_wt.csv")
anbc2.ns.t2 = read.csv("res_ancombc2_NS_T2_wt.csv")

anbc2.as.t1$Taxa= forcats::fct_reorder(anbc2.as.t1$Taxa, anbc2.as.t1$LFC, .desc = TRUE)
anbc2.as.t2$Taxa= forcats::fct_reorder(anbc2.as.t2$Taxa, anbc2.as.t2$LFC, .desc = TRUE)
anbc2.ns.t1$Taxa= forcats::fct_reorder(anbc2.ns.t1$Taxa, anbc2.ns.t1$LFC, .desc = TRUE)
anbc2.ns.t2$Taxa= forcats::fct_reorder(anbc2.ns.t2$Taxa, anbc2.ns.t2$LFC, .desc = TRUE)

comparisons.order = c("FOvsEO", "FOvsIO", "EOvsFO", "IOvsFO")

#custom_palette = c("burlywood","darkseagreen", "springgreen4", "darkorange4")

daa.as.t1 = ggplot(anbc2.as.t1, aes(x = LFC, y = Taxa)) + 
  geom_point(
    aes(fill = Phylum, shape = comparison, color = Phylum), size = 3.5, stroke = 0.1, alpha = 0.6) +
  scale_fill_manual(values = phyla.colors, drop = FALSE) +
  scale_color_manual(values = phyla.colors) +  # Hide color legend, just for fill
  scale_shape_manual(values = c(21, 24)) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        panel.background = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(face = "italic")) +
  ylab("Taxa") +
  xlab("LFC") +
  ggpubr::theme_pubclean() +
  labs(title = "AS T1") +
  geom_vline(xintercept = c(-3, -0.5, 0, 0.5, 3), 
             linetype = c("solid", "dashed", "solid", "dashed", "solid"), 
             size = 0.5) +
  xlim(-3, 3) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  annotate("text", x = -1.5, y = Inf, label = "Forest", vjust = 1, size = 3) +
  annotate("text", x = 1.5, y = Inf, label = "Extensive Intensive", vjust = 1, size = 3)
  
daa.as.t1

daa.as.t2 = ggplot(anbc2.as.t2, aes(x = LFC, y = Taxa)) + 
  geom_point(
    aes(fill = Phylum, shape = comparison, color = Phylum), size = 3.5, stroke = 0.1, alpha = 0.6) +
  scale_fill_manual(values = phyla.colors, drop = FALSE) +
  scale_color_manual(values = phyla.colors) +  # Hide color legend, just for fill
  scale_shape_manual(values = c(21, 24)) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        panel.background = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(face = "italic")) +
  ylab("Taxa") +
  xlab("LFC") +
  ggpubr::theme_pubclean() +
  labs(title = "AS T2") +
  geom_vline(xintercept = 0.5, size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 0.5) + 
  xlim(-4, 4) + 
  geom_vline(xintercept = -4.5, size = 0.5) + 
  geom_vline(xintercept = 3.5, size = 0.5) +
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)) +
  annotate("text", x = -2, y = Inf, label = "Forest", vjust = 1, size = 3) +
  annotate("text", x = 2, y = Inf, label = "Extensive Intensive", vjust = 1, size = 3)

daa.as.t2

both.as = ggpubr::ggarrange(daa.as.t1, daa.as.t2, ncol = 2, labels = c("A"), align = "hv", common.legend = T)
both.as

daa.ns.t1 = ggplot(anbc2.ns.t1, aes(x = LFC, y = Taxa)) + 
  geom_point(
    aes(fill = Phylum, shape = comparison, color = Phylum), size = 3.5, stroke = 0.1, alpha = 0.6) +
  scale_fill_manual(values = phyla.colors, drop = FALSE) +
  scale_color_manual(values = phyla.colors) +  # Hide color legend, just for fill
  scale_shape_manual(values = c(21, 24)) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        panel.background = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(face = "italic")) +
  ylab("Taxa") +
  xlab("LFC") +
  ggpubr::theme_pubclean() +
  labs(title = "NS T1") +
  geom_vline(xintercept = 0.5, size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 0.5) + 
  xlim(-3, 3.5) + 
  geom_vline(xintercept = -3, size = 0.5) + 
  geom_vline(xintercept = 3.5, size = 0.5) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  annotate("text", x = -1.5, y = Inf, label = "Forest", vjust = 1, size = 3) +
  annotate("text", x = 1.75, y = Inf, label = "Extensive Intensive", vjust = 1, size = 3)

daa.ns.t1

daa.ns.t2 = ggplot(anbc2.ns.t2, aes(x = LFC, y = Taxa)) + 
  geom_point(
    aes(fill = Phylum, shape = comparison, color = Phylum), size = 3.5, stroke = 0.1, alpha = 0.6) +
  scale_fill_manual(values = phyla.colors, drop = FALSE) +
  scale_color_manual(values = phyla.colors) +  # Hide color legend, just for fill
  scale_shape_manual(values = c(21, 24)) +
  theme(panel.grid.major = element_line(linetype = 1, color = "grey"),
        panel.background = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(face = "italic")) +
  ylab("Taxa") +
  xlab("LFC") +
  ggpubr::theme_pubclean() +
  labs(title = "NS T2") +
  geom_vline(xintercept = 0.5, size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.5, size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 0.5) + 
  xlim(-2, 4) + 
  geom_vline(xintercept = -2, size = 0.5) + 
  geom_vline(xintercept = 4, size = 0.5) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4)) +
  annotate("text", x = -1, y = Inf, label = "Forest", vjust = 1, size = 3) +
  annotate("text", x = 2, y = Inf, label = "Extensive Intensive", vjust = 1, size = 3)

daa.ns.t2

daa.ns.t1 = daa.ns.t1 + 
  scale_y_discrete(expand = c(0.8, 0.8))
daa.ns.t1

daa.ns.t2 = daa.ns.t2 + 
  scale_y_discrete(expand = c(0.8, 0.8))
daa.ns.t2
both.ns = ggpubr::ggarrange(daa.ns.t1, daa.ns.t2, ncol = 2, labels = c("B"), align = "v", common.legend = T)
both.ns

ggpubr::ggarrange(
  daa.as.t1, daa.as.t2, daa.ns.t1, daa.ns.t2, 
  nrow = 2, ncol = 2, 
  labels = c("A", "B", "C", "D"), 
  align = "h", 
  common.legend = FALSE
)

#Co-occurrence network analyses
#Subsetting the sample
library(NetCoMi)
library(phyloseq)
library(dplyr)
load("phyloseq_obj_biofilms_exp_norm_css.RData")
blik.ott = subset_samples(phy.blik2.css, stream == "otterbach")
blik.pb = subset_samples(phy.blik2.css, stream == "perlenbach")

as.ott = subset_samples(blik.ott, sample_type == "AS")
ns.ott = subset_samples(blik.ott, sample_type == "NS")

as.pb = subset_samples(blik.pb, sample_type == "AS")
ns.pb = subset_samples(blik.pb, sample_type == "NS")

as.ott.in = as.ott %>% subset_samples(land_use == "intensive")
as.ott.in = prune_taxa(taxa_sums(as.ott.in) > 0, as.ott.in)

as.ott.ex = as.ott %>% subset_samples(land_use == "extensive")
as.ott.ex = prune_taxa(taxa_sums(as.ott.ex) > 0, as.ott.ex)

as.ott.fo = as.ott %>% subset_samples(land_use == "forest")
as.ott.fo = prune_taxa(taxa_sums(as.ott.fo) > 0, as.ott.fo)

ns.io = ns.ott %>% subset_samples(land_use == "intensive")
ns.io = prune_taxa(taxa_sums(ns.io) > 0, ns.io)

ns.eo = ns.ott %>% subset_samples(land_use == "extensive")
ns.eo = prune_taxa(taxa_sums(ns.eo) > 0, ns.eo)

ns.fo = ns.ott %>% subset_samples(land_use == "forest")
ns.fo = prune_taxa(taxa_sums(ns.fo) > 0, ns.fo)

as.pb.ex = as.pb %>% subset_samples(land_use == "extensive")
as.pb.ex = prune_taxa(taxa_sums(as.pb.ex) > 0, as.pb.ex)

as.pb.fp = as.pb %>% subset_samples(land_use == "forest")
as.pb.fp = prune_taxa(taxa_sums(as.pb.fp) > 0, as.pb.fp)

ns.pb.ex = ns.pb %>% subset_samples(land_use == "extensive")
ns.pb.ex = prune_taxa(taxa_sums(ns.pb.ex) > 0, ns.pb.ex)

ns.pb.fp = ns.pb %>% subset_samples(land_use == "forest")
ns.pb.fp = prune_taxa(taxa_sums(ns.pb.fp) > 0, ns.pb.fp)


prev_threshold = 0.20  # 20%
nsamples_total = nsamples(as.ott.in)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(as.ott.in), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
as.ott.in.filt = prune_taxa(prev_df >= prev_threshold, as.ott.in)

as.ott.in.glom = tax_glom(as.ott.in.filt, taxrank = "Genus")

as.ott.in.ren = renameTaxa(as.ott.in.glom, 
                            pat = "<name>", 
                            substPat = "<name>_<subst_name>(<subst_R>)",
                            numDupli = "Genus")

net.as.ott.io = netConstruct(as.ott.in.ren,
                         measure = "spieceasi",taxRank = "Genus",
                         verbose = 3, seed = 12345)

nodes = unique(c(net.as.ott.io$edgelist1$v1, net.as.ott.io$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.db.io = netAnalyze(net.as.ott.io, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, 
                           normDeg = FALSE)

sum.net.as.io = summary(props.net.as.io, numbNodes = 5L)
sum.net.as.io$glob_probs
write.csv(sum.net.as2$glob_probs, "net_as_extensive.csv")


conet.db.io = plot(props.net.db.io, 
                 labelScale = FALSE,
                 nodeColor = "cluster", 
                 nodeSize = "eigenvector",
                 nodeFilter = "clustMin",
                 highlightHubs = TRUE,
                 nodeFilterPar = 3,
                 title1 = "AS IO",
                 showTitle = TRUE,
                 cexLabels = 0.5,
                 labelFont = 0.5)

# The taxa you want to keep visible
keep_labels = c("Illumatobacteraceae", "CL500-29 marine group", "Leptothrix", "Rhodoferax", 
                 "Haliangium", "Flectobacillus", "Limnohabitans", "Gemmatimonas", 
                 "Dechloromonas", "Lacihabitans", "Arcicella", 
                 "Candidatus Amoebophilus")

# Grab the original labels from the qgraph object
all_labels = conet.db.io$q1$graphAttributes$Nodes$names

# Create new labels: keep only those in keep_labels, blank the rest
new_labels = ifelse(all_labels %in% keep_labels, all_labels, "")

# Overwrite the labels directly inside the qgraph object
conet.db.io$q1$graphAttributes$Nodes$labels = new_labels

# Now plot normally
pdf("conet_spieceasi_db_io_filt_taxa.pdf", width = 6.41, height = 5.71)

plot(conet.db.io$q1,
     title1 = "DB FO",
     showTitle = TRUE,
     cexLabels = 0.7,
     labelFont = 0.5)

dev.off()

#Co Net of DB EO
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(as.ott.ex)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(as.ott.ex), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
as.ott.ex.filt = prune_taxa(prev_df >= prev_threshold, as.ott.ex)

as.ott.ex.glom = tax_glom(as.ott.ex.filt, taxrank = "Genus")

as.ott.ex.ren = renameTaxa(as.ott.ex.glom, 
                           pat = "<name>", 
                           substPat = "<name>_<subst_name>(<subst_R>)",
                           numDupli = "Genus")

net.as.eo = netConstruct(as.ott.ex.ren,
                             measure = "spieceasi",taxRank = "Genus",
                             verbose = 3, seed = 12345)

nodes = unique(c(net.as.eo$edgelist1$v1, net.as.eo$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.as.eo = netAnalyze(net.as.eo, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)

sum.net.as.eo = summary(props.net.as.eo, numbNodes = 5L)
sum.net.as.eo$glob_probs

conet.db.eo = plot(props.net.as.eo, 
                   labelScale = FALSE,
                   nodeColor = "cluster", 
                   nodeSize = "eigenvector",
                   nodeFilter = "clustMin",
                   highlightHubs = TRUE,
                   nodeFilterPar = 3,
                   title1 = "AS EO",
                   showTitle = TRUE,
                   cexLabels = 0.5,
                   labelFont = 0.5)

# The taxa you want to keep visible
keep_labels = c("Parvibium", "Illumatobacter", "Flavobacterium", "CL500-29 marine group", "Leptothrix",
                "Rhodoferax", "Haliangium", "Flectobacillus", "Gemmatimonas", 
                "Gemmatimonas", "Arenimonas", 
                "Candidatus Amoebophilus")

# Grab the original labels from the qgraph object
all_labels = conet.db.eo$q1$graphAttributes$Nodes$names

# Create new labels: keep only those in keep_labels, blank the rest
new_labels = ifelse(all_labels %in% keep_labels, all_labels, "")

# Overwrite the labels directly inside the qgraph object
conet.db.eo$q1$graphAttributes$Nodes$labels = new_labels

# Now plot normally
pdf("conet_spieceasi_db_eo_filt_taxa.pdf", width = 6.41, height = 5.71)

plot(conet.db.eo$q1,
     title1 = "DB EO",
     showTitle = TRUE,
     cexLabels = 0.7,
     labelFont = 0.5)

dev.off()

#Network of DB FO
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(as.ott.fo)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(as.ott.fo), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
as.ott.fo.filt = prune_taxa(prev_df >= prev_threshold, as.ott.fo)

as.ott.fo.glom = tax_glom(as.ott.fo.filt, taxrank = "Genus")

as.ott.fo.ren = renameTaxa(as.ott.fo.glom, 
                              pat = "<name>", 
                              substPat = "<name>_<subst_name>(<subst_R>)",
                              numDupli = "Genus")

net.db.fo = netConstruct(as.ott.fo.ren,
                       measure = "spieceasi",taxRank = "Genus",
                       verbose = 3, seed = 12345)

#nodes = unique(c(net.db.fo$edgelist1$v1, net.db.fo$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.as.fo = netAnalyze(net.db.fo, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, 
                           normDeg = FALSE)
sum.net.as3 = summary(props.net.as.fo, numbNodes = 5L)
sum.net.as3$hubs
#write.csv(sum.net.as3$glob_probs, "net_as_intensive.csv")

conet.db.fo = plot(props.net.as.fo, 
                 labelScale = FALSE,
                 nodeColor = "cluster", 
                 nodeSize = "eigenvector",
                 nodeFilter = "clustMin",
                 highlightHubs = TRUE,
                 nodeFilterPar = 3,
                 title1 = "DB FO",
                 showTitle = TRUE,
                 cioLabels = 0.5,
                 labelFont = 0.3)

# The taxa you want to keep visible
keep_labels <- c("Sphaerotilus", "Thermomonas", "Illumatobacter", "Parvibium", 
                 "Deinococcaceae", "Acidimicrobiia", "Candidatus Amoebophilus", 
                 "CL500-29 marine group", "Gemmatimonas", "Leptothrix", 
                 "Rhodoferax", "Haliangium")

# Grab the original labels from the qgraph object
all_labels <- conet.db.fo$q1$graphAttributes$Nodes$names

# Create new labels: keep only those in keep_labels, blank the rest
new_labels <- ifelse(all_labels %in% keep_labels, all_labels, "")

# Overwrite the labels directly inside the qgraph object
conet.db.fo$q1$graphAttributes$Nodes$labels <- new_labels

# Now plot normally
pdf("filtered_network_plot.pdf", width = 6.41, height = 5.71)

plot(conet.db.fo$q1,
     title1 = "DB FO",
     showTitle = TRUE,
     cexLabels = 0.7,
     labelFont = 0.5)

dev.off()
#Network of MB FO
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(ns.fo)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(ns.fo), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
ns.ott.fo.filt = prune_taxa(prev_df >= prev_threshold, ns.fo)

ns.ott.fo.glom = tax_glom(ns.ott.fo.filt, taxrank = "Genus")

ns.ott.fo.ren = renameTaxa(ns.ott.fo.glom, 
                           pat = "<name>", 
                           substPat = "<name>_<subst_name>(<subst_R>)",
                           numDupli = "Genus")

net.ns.fo = netConstruct(ns.ott.fo.ren,
                         measure = "spieceasi", taxRank = "Genus",
                         sparsMethod = "t-test",
                         verbose = 3, seed = 12345)

#nodes = unique(c(net.ns.fo$edgelist1$v1, net.ns.fo$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.ns.fo = netAnalyze(net.ns.fo, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)

sum.net.ns.fo = summary(props.net.ns.fo, numbNodes = 5L)
sum.net.ns.fo$hubs
write.csv(sum.net.ns.fo$glob_probs, "net_ns_forest.csv")

conet.mb.fo = plot(props.net.ns.fo, 
                   labelScale = FALSE,
                   nodeColor = "cluster", 
                   nodeSize = "eigenvector",
                   nodeFilter = "clustMin",
                   highlightHubs = TRUE,
                   nodeFilterPar = 3,
                   title1 = "MB FO",
                   showTitle = TRUE,
                   cioLabels = 0.5,
                   labelFont = 0.3)

library(qgraph)

keep_labels = c("Nitrospira", "Lysobacter", "Cloacibacterium", "Ellin6067")

# Grab the original labels from the qgraph object
all_labels = conet.mb.fo$q1$graphAttributes$Nodes$names

# Create new labels: keep only those in keep_labels, blank the rest
new_labels <- ifelse(all_labels %in% keep_labels, all_labels, "")

# Overwrite the labels directly inside the qgraph object
conet.mb.fo$q1$graphAttributes$Nodes$labels <- new_labels

# Now plot normally
pdf("conet_spieceasi_mb_fo_filt_taxa.pdf", width = 6.41, height = 5.71)

plot(conet.mb.fo$q1,
     title1 = "MB FO",
     showTitle = TRUE,
     cexLabels = 0.7,
     labelFont = 0.5)
dev.off()

#Network IO NS
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(ns.io)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(ns.io), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
ns.io.filt = prune_taxa(prev_df >= prev_threshold, ns.io)

ns.io.filt = tax_glom(ns.io.filt, taxrank = "Genus")

ns.io.filt.ren = renameTaxa(ns.io.filt, 
                                pat = "<name>", 
                                substPat = "<name>_<subst_name>(<subst_R>)",
                                numDupli = "Genus")

net.ns.io = netConstruct(ns.io.filt.ren,
                       measure = "spieceasi",taxRank = "Genus",
                       sparsMethod = "t-test",
                       verbose = 3, seed = 12345)

nodes = unique(c(net.ns.io$edgelist1$v1, net.ns.io$edgelist1$v2)) %>% as.data.frame()

props.net.ns.io = netAnalyze(net.ns.io, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)

sum.net.net.ns.io = summary(props.net.ns.io, numbNodes = 5L)

sum.net.net.ns.io$hubs
write.csv(sum.net.net.ns.io$glob_probs, "net_ns_intensive.csv")
par(cex = 0.5)
net.ns.io = plot(props.net.ns.io, 
                 labelScale = FALSE,
                 nodeColor = "cluster", 
                 nodeSize = "eigenvector",
                 nodeFilter = "clustMin",
                 highlightHubs = TRUE,
                 nodeFilterPar = 3,
                 title1 = "MB IO",
                 showTitle = TRUE,
                 cioLabels = 0.5,
                 labelFont = 0.5)

# The taxa you want to keep visible
keep_labels = c("Cloacibacterium", "Ferruginibacter", "Blastocatella", "Nitrospira")

# Grab the original labels from the qgraph object
all_labels = net.ns.io$q1$graphAttributes$Nodes$names

# Create new labels: keep only those in keep_labels, blank the rest
new_labels = ifelse(all_labels %in% keep_labels, all_labels, "")

# Overwrite the labels directly inside the qgraph object
net.ns.io$q1$graphAttributes$Nodes$labels <- new_labels

# Now plot normally
pdf("filtered_network_plot.pdf", width = 6.41, height = 5.71)

plot(net.ns.io$q1,
     title1 = "MB IO",
     showTitle = TRUE,
     cexLabels = 0.7,
     labelFont = 0.5)

dev.off()

#Network EO NS
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(ns.eo)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(ns.eo), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
ns.eo.filt = prune_taxa(prev_df >= prev_threshold, ns.eo)

ns.eo.filt = tax_glom(ns.eo.filt, taxrank = "Genus")

ns.eo.filt.ren = renameTaxa(ns.eo.filt, 
                            pat = "<name>", 
                            substPat = "<name>_<subst_name>(<subst_R>)",
                            numDupli = "Genus")

net.ns.eo = netConstruct(ns.eo.filt.ren,
                         measure = "spieceasi",taxRank = "Genus",
                         verbose = 3, seed = 12345)

props.net.ns.eo = netAnalyze(net.ns.eo, 
                            centrLCC = TRUE,
                            clustMethod = "cluster_fast_greedy",
                            hubPar = "eigenvector",
                            weightDeg = FALSE, 
                            normDeg = FALSE)

sum.net.ns.eo = summary(props.net.ns.eo, numbNodes = 5L)

sum.net.ns.eo$glob_probs
write.csv(sum.net.ns3$glob_probs, "net_ns_extensive.csv")

nodes = unique(c(net.ns.eo$edgelist1$v1, net.ns.eo$edgelist1$v2)) %>% as.data.frame()

par(cex = 0.7)  # globally reduce all text
net.ns.eo = plot(props.net.ns.eo, 
                 labelScale = FALSE,
                 nodeColor = "cluster", 
                 nodeSize = "eigenvector",
                 nodeFilter = "clustMin",
                 highlightHubs = TRUE,
                 nodeFilterPar = 3,
                 title1 = "NS EO",
                 showTitle = TRUE,
                 labelSize = 0.3,
                 labelFont = 0)

#Perlenbach networks

#Network of DB FP
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(as.pb.fp)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(as.pb.fp), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
as.pb.fp.filt = prune_taxa(prev_df >= prev_threshold, as.pb.fp)

as.pb.fp.glom = tax_glom(as.pb.fp.filt, taxrank = "Genus")

as.pb.fp.ren = renameTaxa(as.pb.fp.glom, 
                           pat = "<name>", 
                           substPat = "<name>_<subst_name>(<subst_R>)",
                           numDupli = "Genus")

net.db.fp = netConstruct(as.pb.fp.ren,
                         measure = "spieceasi",taxRank = "Genus",
                         verbose = 3, seed = 12345)

nodes = unique(c(net.db.fp$edgelist1$v1, net.db.fp$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.db.fp = netAnalyze(net.db.fp, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)
sum.net.db.fp = summary(props.net.db.fp, numbNodes = 5L)
sum.net.db.fp$glob_probs
sum.net.db.fp$hubs
#write.csv(sum.net.as3$glob_probs, "net_as_intensive.csv")

conet.db.fo = plot(props.net.db.fp, 
                   labelScale = FALSE,
                   labels = FALSE,
                   nodeColor = "cluster", 
                   nodeSize = "eigenvector",
                   nodeFilter = "clustMin",
                   highlightHubs = TRUE,
                   nodeFilterPar = 3,
                   title1 = "DB FP",
                   showTitle = TRUE,
                   cioLabels = 0.5,
                   labelFont = 0.3)

#Network of DB EP
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(as.pb.ex)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(as.pb.ex), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
as.pb.ex.filt = prune_taxa(prev_df >= prev_threshold, as.pb.ex)

as.pb.ex.glom = tax_glom(as.pb.ex.filt, taxrank = "Genus")

as.pb.ex.ren = renameTaxa(as.pb.ex.glom, 
                          pat = "<name>", 
                          substPat = "<name>_<subst_name>(<subst_R>)",
                          numDupli = "Genus")

net.db.ep = netConstruct(as.pb.ex.ren,
                         measure = "spieceasi",taxRank = "Genus",
                         verbose = 3, seed = 12345)

nodes = unique(c(net.db.ep$edgelist1$v1, net.db.ep$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.db.ep = netAnalyze(net.db.ep, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)
props.net.db.ep = summary(props.net.db.ep, numbNodes = 5L)
props.net.db.ep$glob_probs
props.net.db.ep$hubs
#write.csv(sum.net.as3$glob_probs, "net_as_intensive.csv")

conet.db.ep = plot(props.net.db.ep, 
                   labelScale = FALSE,
                   labels = FALSE,
                   nodeColor = "cluster", 
                   nodeSize = "eigenvector",
                   nodeFilter = "clustMin",
                   highlightHubs = TRUE,
                   nodeFilterPar = 3,
                   title1 = "DB EP",
                   showTitle = TRUE,
                   cioLabels = 0.5,
                   labelFont = 0.3)


#CONETs of MB PB samples
#Network of MB EP
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(ns.pb.ex)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(ns.pb.ex), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
ns.pb.ex.filt = prune_taxa(prev_df >= prev_threshold, ns.pb.ex)

ns.pb.ex.glom = tax_glom(ns.pb.ex.filt, taxrank = "Genus")

ns.pb.ex.ren = renameTaxa(ns.pb.ex.glom, 
                          pat = "<name>", 
                          substPat = "<name>_<subst_name>(<subst_R>)",
                          numDupli = "Genus")

net.mb.ep = netConstruct(ns.pb.ex.ren,
                         measure = "spieceasi", taxRank = "Genus",
                         verbose = 3, seed = 12345)

nodes = unique(c(net.mb.ep$edgelist1$v1, net.mb.ep$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.mb.ep = netAnalyze(net.mb.ep, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)

sum.net.mb.ep = summary(props.net.mb.ep, numbNodes = 5L)
sum.net.mb.ep$glob_probs
sum.net.mb.ep$hubs
#write.csv(sum.net.ns3$glob_probs, "net_ns_intensive.csv")

conet.mb.ep = plot(props.net.mb.ep, 
                   labelScale = FALSE,
                   labels = FALSE,
                   nodeColor = "cluster", 
                   nodeSize = "eigenvector",
                   nodeFilter = "clustMin",
                   highlightHubs = TRUE,
                   nodeFilterPar = 3,
                   title1 = "MB EP",
                   showTitle = TRUE,
                   cioLabels = 0.5,
                   labelFont = 0.3)

#Network of MB FP
prev_threshold = 0.20  # 20%
nsamples_total = nsamples(ns.pb.fp)

# Calculate prevalence of each taxon
prev_df = apply(otu_table(ns.pb.fp), 1, function(x) sum(x > 0) / nsamples_total)

# Keep only taxa with ≥20% prevalence
ns.pb.fp.filt = prune_taxa(prev_df >= prev_threshold, ns.pb.fp)

ns.pb.fp.glom = tax_glom(ns.pb.fp.filt, taxrank = "Genus")

ns.pb.fp.ren = renameTaxa(ns.pb.fp.glom, 
                          pat = "<name>", 
                          substPat = "<name>_<subst_name>(<subst_R>)",
                          numDupli = "Genus")

net.mb.fp = netConstruct(ns.pb.fp.ren,
                         measure = "spieceasi", taxRank = "Genus",
                         verbose = 3, seed = 12345)

nodes = unique(c(net.mb.fp$edgelist1$v1, net.mb.fp$edgelist1$v2)) %>% as.data.frame() #to get the number of nodes

props.net.mb.fp = netAnalyze(net.mb.fp, 
                             centrLCC = TRUE,
                             clustMethod = "cluster_fast_greedy",
                             hubPar = "eigenvector",
                             weightDeg = FALSE, 
                             normDeg = FALSE)

sum.props.net.mb.fp = summary(props.net.mb.fp, numbNodes = 5L)
sum.props.net.mb.fp$glob_probs
sum.props.net.mb.fp$hubs
#write.csv(sum.net.ns3$glob_probs, "net_ns_intensive.csv")

conet.mb.fp = plot(props.net.mb.fp, 
                   labelScale = FALSE,
                   labels = FALSE,
                   nodeColor = "cluster", 
                   nodeSize = "eigenvector",
                   nodeFilter = "clustMin",
                   highlightHubs = TRUE,
                   nodeFilterPar = 3,
                   title1 = "MB FP",
                   showTitle = TRUE,
                   cioLabels = 0.5,
                   labelFont = 0.3)

###rain plot of Fig. S3
rain = read.csv("niederschlag_2023_eiserzell.csv")
library(ggplot2)

rain$Datum = as.Date(rain$Datum, format = "%d/%m/%Y")
ggplot() + geom_bar(rain, mapping = aes(x = Datum, y = avg), fill = "grey", alpha = 0.4,
                    stat = "identity") +
  geom_line(rain, mapping = aes(x = Datum, y = Niederschlag), col = "royalblue4") +
  scale_x_date(
    date_labels = "%b",  # Show only month names
    breaks = "1 month"   # Place ticks at monthly intervals
  ) +
  scale_y_continuous(name = "") + ylim(0, 200) + ggpubr::theme_pubr() + labs(x = "2023", y = "Daily average and accumulated rainfall (mm)")


###DOC plot of Fig. S4
doc.df = read.csv("DOC_data_sebastian.csv", header = T)
# Clean timestamps (handle case with only date)
doc.df = doc.df %>%
  mutate(Timestamp = as.character(Timestamp)) %>% # ensure it's character
  separate(Timestamp, into = c("Date", "Hour"), sep = " ", fill = "right")

doc.df[is.na(doc.df)] <- 0

# Optional: pivot longer to have one variable per line (IO and FO)
df.long = pivot_longer(doc.df, cols = c("IO", "FO"), names_to = "Type", values_to = "Value")
daily_avg = doc.long %>%
  group_by(Date, Month, Type) %>%
  summarise(Daily_Avg = mean(Value, na.rm = TRUE)) %>%
  ungroup()

daily_avg$Date = as.Date(daily_avg$Date, format = "%d/%m/%Y")  # Adjust format as needed

# Create the plot
scat.doc = ggplot(daily_avg, aes(x = Date, y = Daily_Avg, color = Type, group = Type)) +
  geom_point(size = 1) +
  geom_line(size = 1) +
  scale_color_manual(values = c("FO" = "#009E73", "IO" = "#E69F00")) +
  scale_x_date(
    date_labels = "%b",  # Show only month names
    breaks = "1 month"   # Place ticks at monthly intervals
  ) +
  labs(
    x = "Month",
    y = "Average DOC concentration [mg/L]",
    color = "Site"
  ) + ylim(0,40) +
  ggpubr::theme_pubr()

scat.doc

conet.as.ott.in.t1 = SpiecEasi::spiec.easi(as.ott.in.t1)


######PCoA with Taxa as vectors
library(vegan)
phy.genus = tax_glom(as.ott, taxrank = "Genus")
bray.as.ott = phyloseq::distance(as.ott, method = "bray")
pcoa.res = ape::pcoa(bray.as.ott)
pcoa.as.ott.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                            Axis1 = pcoa.res$vectors[, 1],
                            Axis2 = pcoa.res$vectors[, 2])

sd.ott.as = read.csv("sd_as_ott.csv", row.names = 1)


pcoa.df = merge(pcoa.as.ott.df, sd.ott.as, by.x = "SampleID", by.y = "row.names")
pcoa.df$SampleID.y = NULL
pcoa.df$sample_type = NULL
#pcoa.df$SampleID = NULL

phy.genus = tax_glom(as.ott, taxrank = "genus")
asvs.as.ott = otu_table(as.ott) %>% as.data.frame() %>% t()
all(rownames(sample_data(as.ott)) %in% rownames(asvs.as.ott))

asvs.hel.as.ott = decostand(asvs.as.ott, method = "hellinger")
asvs.var = apply(asvs.hel.as.ott, 2, var)
asvs.hel.filtered = asvs.hel.as.ott[, order(asvs.var, decreasing = TRUE)[1:5000]]
env.fit = envfit(pcoa.res$vectors, asvs.hel.filtered, permutations = 999)
# Extract p-values
pvals = env.fit$vectors$pvals
padj = p.adjust(pvals, method = "fdr")
signif.vars = names(padj[padj < 0.05])

# Filter the vectors
significant.vectors = env.fit$vectors$arrows[signif.vars, ]

env.vectors = as.data.frame(scores(significant.vectors))
env.vectors$Variable = rownames(env.vectors)

lu.colors = c("intensive" = "#E69F00", "forest" = "#009E73", "extensive" = "#999999")
land.use = sample_data(as.ott)$land_use
st = sample_data(as.ott)$st
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Compute length only for the significant vectors
env.vectors$Length <- env.fit$vectors$r[rownames(env.vectors)]
env.vectors$x_end <- env.vectors$Axis.1 * env.vectors$Length * 0.5
env.vectors$y_end <- env.vectors$Axis.2 * env.vectors$Length * 0.5

# Get Genus from tax_table only for significant taxa
tax_table = as.data.frame(tax_table(as.ott))
env.vectors$Genus <- tax_table[rownames(env.vectors), "Genus"]

# Filter out NA Genus
env.vectors <- env.vectors[!is.na(env.vectors$Genus), ]
env.vectors <- env.vectors[env.vectors$Length > 0.75, ]

tax_table = as.data.frame(tax_table(as.ott))
env.vectors$Genus <- tax_table[rownames(env.vectors), "Genus"]
env.vectors <- env.vectors[!is.na(env.vectors$Genus), ]
# Then plot:
pcoa.as.ott = ggplot(pcoa.df, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(land_use), shape = factor(st)), size = 3, alpha = 0.6) +
  scale_color_manual(values = lu.colors) +
  geom_segment(data = env.vectors,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = env.vectors,
                  aes(x = x_end, y = y_end, label = Genus),
                  size = 4) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA AS") +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  theme_bw()
pcoa.as.ott

###PCoA PB AS with Taxa
blik.pb = subset_samples(phy.blik2.css, stream == "perlenbach")
as.pb = subset_samples(blik.pb, sample_type == "AS")

bray.as.pb = phyloseq::distance(as.pb, method = "bray")
pcoa.res = ape::pcoa(bray.as.pb)
pcoa.as.pb.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                            Axis1 = pcoa.res$vectors[, 1],
                            Axis2 = pcoa.res$vectors[, 2])

sd.pb.as = read.csv("sd_pb_as.csv", row.names = 1)


pcoa.df = merge(pcoa.as.pb.df, sd.pb.as, by.x = "SampleID", by.y = "row.names")
pcoa.df$SampleID.y = NULL
pcoa.df$sample_type = NULL
#pcoa.df$SampleID = NULL

asvs.as.pb = otu_table(as.pb) %>% as.data.frame() %>% t()
all(rownames(sample_data(as.pb)) %in% rownames(asvs.as.pb))

asvs.hel.as.pb = decostand(asvs.as.pb, method = "hellinger")
asvs.var = apply(asvs.hel.as.pb, 2, var)
asvs.hel.filtered = asvs.hel.as.pb[, order(asvs.var, decreasing = TRUE)[1:5000]]
env.fit = envfit(pcoa.res$vectors, asvs.hel.filtered, permutations = 999)
# Extract p-values
pvals = env.fit$vectors$pvals
padj = p.adjust(pvals, method = "fdr")
signif.vars = names(padj[padj < 0.05])

# Filter the vectors
significant.vectors = env.fit$vectors$arrows[signif.vars, ]

env.vectors = as.data.frame(scores(significant.vectors))
env.vectors$Variable = rownames(env.vectors)

lu.colors = c("intensive" = "#E69F00", "forest" = "#009E73", "extensive" = "#999999")
land.use = sample_data(as.ott)$land_use
st = sample_data(as.ott)$st
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Compute length only for the significant vectors
env.vectors$Length <- env.fit$vectors$r[rownames(env.vectors)]
env.vectors$x_end <- env.vectors$Axis.1 * env.vectors$Length * 0.5
env.vectors$y_end <- env.vectors$Axis.2 * env.vectors$Length * 0.5

# Get Genus from tax_table only for significant taxa
tax_table = as.data.frame(tax_table(as.pb))
env.vectors$Genus <- tax_table[rownames(env.vectors), "Genus"]

# Filter out NA Genus
env.vectors <- env.vectors[!is.na(env.vectors$Genus), ]
env.vectors <- env.vectors[env.vectors$Length > 0.75, ]

# Then plot:
pcoa.as.pb = ggplot(pcoa.df, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(land_use), shape = factor(st)), size = 3, alpha = 0.6) +
  scale_color_manual(values = lu.colors) +
  geom_segment(data = env.vectors,
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = env.vectors,
                  aes(x = x_end, y = y_end, label = Genus),
                  size = 4) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA AS") +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, by = 0.25)) +
  theme_bw()
pcoa.as.pb
