---
  title: "iCAMP Biofilms experiment"
output: html_document
source: "https://github.com/DaliangNing/iCAMP1"
---
  
library(iCAMP)
library(phyloseq)
library(dplyr)
library(ape)
library(ggplot2)
library(wesanderson)
library(ggpubr)
library(ggpubr)
library(grid)

#Load phyloseq object and subset
blik.ott = subset_samples(phy.blik2.css, river == "otterbach")
as.ott = subset_samples(blik.ott, sample_type == "DB")
ns.ott = subset_samples(blik.ott, sample_type == "MB")

as.ott.t1 = subset_samples(as.ott, st == "T1")
as.ott.t2 = subset_samples(as.ott, st == "T2")

as.ott.t1 = prune_taxa(taxa_sums(as.ott.t1) > 0, as.ott.t1)
as.ott.t2 = prune_taxa(taxa_sums(as.ott.t2) > 0, as.ott.t2)

ns.ott.t1 = subset_samples(ns.ott, st == "T1")
ns.ott.t2 = subset_samples(ns.ott, st == "T2")
ns.ott.t1 = prune_taxa(taxa_sums(ns.ott.t1) > 0, ns.ott.t1)
ns.ott.t2 = prune_taxa(taxa_sums(ns.ott.t2) > 0, ns.ott.t2)

blik.pb = subset_samples(phy.blik2.css, stream == "perlenbach")
as.pb = subset_samples(blik.pb, sample_type == "DB")
ns.pb = subset_samples(blik.pb, sample_type == "MB")

as.pb.t1 = subset_samples(as.pb, st == "T1")
as.pb.t2 = subset_samples(as.pb, st == "T2")
as.pb.t1 = prune_taxa(taxa_sums(as.pb.t1) > 0, as.pb.t1)
as.pb.t2 = prune_taxa(taxa_sums(as.pb.t2) > 0, as.pb.t2)

ns.pb.t1 = subset_samples(ns.pb, st == "T1")
ns.pb.t2 = subset_samples(ns.pb, st == "T2")
ns.pb.t1 = prune_taxa(taxa_sums(ns.pb.t1) > 0, ns.pb.t1)
ns.pb.t2 = prune_taxa(taxa_sums(ns.pb.t2) > 0, ns.pb.t2)

###Follow the these steps with each of the phyloseq subsets
#load required df and shuffle them a bit
otu = otu_table(as.ott.t1) %>% as.data.frame()
write.csv(otu, "com_as_ott_t1.csv")
treat2col = sample_data(ns.pb.t2) %>% as.data.frame()
env = read.csv("environment_as_ott_t1.csv", row.names = 1)

treat2col$SampleID = NULL
treat2col$sample = NULL
treat2col$stream = NULL 
treat2col$sample_type = NULL 
treat2col$site = NULL
treat2col$st = NULL
treat2col$norm_factor = NULL

treat = treat2col

treat$SampleID = row.names(treat)

wd=""
save.wd="AS_OTT_T1"
if(!dir.exists(save.wd)){dir.create(save.wd)}

prefix = "AS_OTT_T1"  # prefix of the output file names. usually use a project ID.
rand.time = 100  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker = 4 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G = 50 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

tree = phy_tree(ns.pb.t2)
write.tree(tree, file = "ns_pb_t2.nwk")
tree.file = "ns_pb_t2.nwk"
tree=read.tree(file = tree.file)

comm=t(read.csv("com_ns_pb_t2.csv", header = T, sep = ",", row.names = 1, as.is = T, stringsAsFactors = F, comment.char = "", check.names = F))

comm = t(otu) %>% as.data.frame()

clas = tax_table(ns.pb.t2) %>% as.data.frame()
clas$Species.1 = NULL

env=read.table("environment_as_ott_t1.csv", header = TRUE, sep = ",", row.names = 1,
               as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
               check.names = FALSE) # skip this if you do not have env.file

treat = as.data.frame(treat)
# 4 # match sample IDs in OTU table and treatment information table
sampid.check=match.name(rn.list=list(comm=comm, treat=treat2col, env=env))
# sampid.check=match.name(rn.list=list(comm=comm,treat=treat)) # if you do not have env.file
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched their IDs, the unmatched samples will be removed.
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
env=sampid.check$env # skip this if you do not have env.file

# 5 # match OTU IDs in OTU table and tree file
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(clas=clas),tree.list=list(tree=tree))
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched the IDs before, the unmatched OTUs will be removed.
comm=spid.check$comm
clas=spid.check$clas
tree=spid.check$tree

# 6 # calculate pairwise phylogenetic distance matrix.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large phylogenetic distance matrix. 
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G =memory.G)
  # output files:
  # path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
  # pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
  # pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
  # pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.
}else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

# 7 # assess niche preference difference between species
# env is required for this step.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large niche difference matrix. 
niche.dif=iCAMP::dniche(env = env, comm = comm, method = "niche.value",
                        nworker = nworker, out.dist=FALSE, bigmemo=TRUE,
                        nd.wd=save.wd)

# 8 # within-bin phylogenetic signal assessment.
# For real data, you may try several different settings of binning, and choose the one leading to the best within-bin phylogenetic signal.
# env is required for this step.
# 8.1 # try phylogenetic binning using current setttings.
ds = 0.2 # setting can be changed to explore the best choice
bin.size.limit = 12 # setting can be changed to explore the best choice. # here set as 5 just for the small example dataset. For real data, usually try 12 to 48.

# The tree for taxa.binphy.big must be a rooted tree.
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)
  tree=tree.rt$tree
}
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)

sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=3 # you may remove some species, if they are too rare to perform reliable correlation test.
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)
if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)

bin.size.limit = 12 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
sig.index="Confidence" # see other options in help document of icamp.big.
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no", meta.ab = NULL)

#bin level statistics
treat2=data.frame(CountAll="All", treat, stringsAsFactors = FALSE)
icbin=icamp.bins(icamp.detail = icres$detail, treat = treat2,
                 clas = clas, boot = TRUE, rand.time = 1000)
save(icbin, file = paste0(prefix,".iCAMPBinSummaryDetail.rda"))

Ptk=icbin$Ptk
BPtk=icbin$BPtk
BRPtk=icbin$BRPtk
clas.bin=icbin$Class.Bin
bin.clas=icbin$Bin.TopClass
Binwt=icbin$Binwt
save(Ptk, prefix=prefix, file = "ProcessImpBin.Ptk")
save(BPtk,prefix=prefix,file = "BinContrProcess.BPtk")
save(BRPtk,prefix=prefix,file = "BinRelContrProcess.BPtk")
save(clas.bin,prefix=prefix,file = "OTUBinClass")
save(bin.clas,prefix=prefix,file = "BinTopOTU")
save(Binwt,prefix=prefix,file = "BinRelAbundance")

treat2 = treat
treat2$SampleID = row.names(treat2)
treat2 = treat2[, c(ncol(treat2), 1:(ncol(treat2) - 1))]

icres.detail = as.matrix(icres$detail)

treat2 = treat2[treat2$SampleID %in% row.names(icres$detail), ]
treat2 = treat2[match(row.names(icres$detail), treat2$SampleID), ]

icbin=iCAMP::icamp.bins(icamp.detail = icres$detail,treat = treat,
                        clas=clas, silent=FALSE, boot = TRUE,
                        rand.time = rand.time, between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
write.csv(icbin$Pt ,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)

caa.as.pb.t1 <- read.csv("AS_PB_T1.ProcessImportance_EachGroup.csv") %>% .[-1, ] %>% mutate(sample = "DB") %>% mutate(st = "T1")
caa.as.pb.t2 = read.csv("AS_PB_T2.ProcessImportance_EachGroup.csv") %>% .[-1, ] %>% mutate(sample = "DB") %>% mutate(st = "T2")
caa.ns.pb.t1 = read.csv("NS_PB_T1.ProcessImportance_EachGroup.csv") %>% .[-1, ] %>% mutate(sample = "MB") %>% mutate(st = "T1")
caa.ns.pb.t2 = read.csv("NS_PB_T2.ProcessImportance_EachGroup.csv") %>% .[-1, ] %>% mutate(sample = "MB") %>% mutate(st = "T2")

caa.pb = rbind(caa.as.pb.t1, caa.as.pb.t2, caa.ns.pb.t1, caa.ns.pb.t2)
caa.pb.melt = reshape2::melt(caa.pb)
caa.pb.melt$percentage = caa.pb.melt$value * 100
caa.pb.melt$land_use = caa.pb.melt$Group

icamp = read.csv("iCAMP_NS_AS.csv")
icamp_melt = reshape2::melt(icamp)
icamp_melt$value100 = icamp_melt$value * 100
icamp_melt$value = NULL
icamp_melt$value100 = icamp_melt$value

icamp.as = subset.data.frame(icamp_melt, Sample_type == "DB")
icamp.ns = subset.data.frame(icamp_melt, Sample_type == "MB")

assembly_palette = wesanderson::wes_palette(name = "Cavalcanti1")
icamp.as$land_use = factor(icamp.as$land_use, levels = c("FO", "EO", "IO"))
icamp.ns$land_use = factor(icamp.ns$land_use, levels = c("FO", "EO", "IO"))

ott.order = c("EO", "IO", "FO")

icamp.t1.sb = ggplot(icamp.as, aes(x = factor(land_use, ordered = ott.order), y = value100, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(value100, 1), "%")), 
            position = position_stack(vjust = 0.5), color = "white") +
  facet_grid(~ st, scales = "free_y") +
  scale_fill_manual(values = assembly_palette) +
  labs(y = "Percentage", x = "AS") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "top")

icamp.t1.sb

icamp.t2.sb = ggplot(icamp.ns, aes(x = land_use, y = value100, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(value100, 1), "%")), 
            position = position_stack(vjust = 0.5), color = "white") +
  facet_grid(~ st, scales = "free_y") +
  scale_fill_manual(values = assembly_palette) +
  labs(y = "Percentage", x = "NS") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "top")

caa.as.ns.ott = ggpubr::ggarrange(icamp.t1.sb, icamp.t2.sb, ncol = 2, labels=c("A"), common.legend = T)
caa.as.ns.ott.annot = annotate_figure(
  caa.as.ns.ott,
  right = text_grob(
    "Otterbach", 
    rot = 270, 
    size = 11, 
    face = "plain", 
    family = "Arial"
  ))
caa.as.ns.ott.annot

caa.as.pb = subset.data.frame(caa.pb.melt, sample == "DB")
caa.ns.pb = subset.data.frame(caa.pb.melt, sample == "MB")

assembly_palette = wesanderson::wes_palette(name = "Cavalcanti1")
caa.as.pb$land_use = factor(caa.as.pb$land_use, levels = c("forest", "extensive"))
caa.ns.pb$land_use = factor(caa.ns.pb$land_use, levels = c("forest", "extensive"))

plot.caa.as.pb = ggplot(caa.as.pb, aes(x = land_use, y = percentage, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), color = "white") +
  facet_grid(~ st, scales = "free_y") +
  scale_fill_manual(values = assembly_palette) +
  labs(y = "Percentage", x = "AS") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "top")
plot.caa.as.pb

plot.caa.ns.pb = ggplot(caa.ns.pb, aes(x = land_use, y = percentage, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), color = "white") +
  facet_grid(~ st, scales = "free_y") +
  scale_fill_manual(values = assembly_palette) +
  labs(y = "Percentage", x = "NS") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "top")
plot.caa.ns.pb

caa.as.ns.pb = ggpubr::ggarrange(
  plot.caa.as.pb, plot.caa.ns.pb, 
  ncol = 2, 
  labels = c("B"), 
  align = "hv", 
  common.legend = TRUE
)
caa.as.ns.pb.annot = annotate_figure(
  caa.as.ns.pb,
  right = text_grob(
    "Perlenbach", 
    rot = 270, 
    size = 11, 
    face = "plain", 
    family = "Arial"
  ))
caa.as.ns.pb.annot

caa.as.ns.ott.ob = ggpubr::ggarrange(
  caa.as.ns.ott.annot, caa.as.ns.pb.annot, 
  nrow = 2, 
  align = "hv", 
  common.legend = TRUE
)
caa.as.ns.ott.ob