#### 0. Load required packages and data ####
# install.packages("BiocManager")
library(phyloseq)
library(vegan)
library(GUniFrac)

get.unifrac <- function(obj.phylo) {
  tr <- phy_tree(obj.phylo)
  tmp.otu <- as.data.frame(t(otu_table(obj.phylo)@.Data))
  tmp.ret <- GUniFrac(tmp.otu, tr)
  cat("finish computing the unifrac distance \n")
  tmp.d = tmp.ret$unifracs[, , "d_UW"]
  metap.nmds = vegan::metaMDS(tmp.d, trymax = 100)
}

set.seed(42)
setwd("./data")
load("taxa_name_NorTA.Rdata", verbose = T)  # Load taxon names
load("phylo_tree_info.Rdata", verbose = T)  # Load phylogenetic tree
load("unifrac_nmds.Rdata")  # Load existing NMDS results (optional)

#### 1. Process CASE group data (order: real → MBDDPM → WGAN → MIDASim → SparseDOSSA2 → NorTA) ####
load("Data_NielsenHB_real.Rdata", verbose = T)  # Real data
load("Data_by_MBDDPM_v1.0.Rdata", verbose = T)  # MBDDPM
load("Data_by_MBDDPM_NONE_v1.0.Rdata", verbose = T)  # MBDDPM
load("Data_by_WGAN_v1.0.Rdata", verbose = T)    # WGAN
load("Data_by_MIDASim_v1.0.Rdata", verbose = T) # MIDASim
load("Data_by_SparseDOSSA2_v1.0.Rdata", verbose = T) # SparseDOSSA2
load("Data_by_NorTA.Rdata")

# (1) MBDDPM
case.MBDDPM.filter = case.MBDDPM[, colnames(case.MBDDPM) %in% taxa.case.NorTA]
case.MBDDPM.tab = otu_table(t(case.MBDDPM.filter), taxa_are_rows = T)
case.MBDDPM.phylo = phyloseq(case.MBDDPM.tab, tree.tip)
nmds.case.MBDDPM = get.unifrac(case.MBDDPM.phylo)

# MBDDPM_NONE
case.MBDDPM_NONE.filter = case.MBDDPM_NONE[, colnames(case.MBDDPM_NONE) %in% taxa.case.NorTA]
case.MBDDPM_NONE.tab = otu_table(t(case.MBDDPM_NONE.filter), taxa_are_rows = T)
case.MBDDPM_NONE.phylo = phyloseq(case.MBDDPM_NONE.tab, tree.tip)
nmds.case.MBDDPM_NONE = get.unifrac(case.MBDDPM_NONE.phylo)

# (2) WGAN
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]
case.WGAN.tab = otu_table(t(case.WGAN.filter), taxa_are_rows = T)
case.WGAN.phylo = phyloseq(case.WGAN.tab, tree.tip)
nmds.case.WGAN = get.unifrac(case.WGAN.phylo)

# (3) MIDASim
case.MIDASim.filter = case.MIDASim[, colnames(case.MIDASim) %in% taxa.case.NorTA]
case.MIDASim.tab = otu_table(t(case.MIDASim.filter), taxa_are_rows = T)
case.MIDASim.phylo = phyloseq(case.MIDASim.tab, tree.tip)
nmds.case.MIDASim = get.unifrac(case.MIDASim.phylo)

# (4) SparseDOSSA2
case.SparseDOSSA.filter = case.compo.SparseDOSSA2[, colnames(case.compo.SparseDOSSA2) %in% taxa.case.NorTA]
case.SparseDOSSA.tab = otu_table(t(case.SparseDOSSA.filter), taxa_are_rows = T)
case.SparseDOSSA.phylo = phyloseq(case.SparseDOSSA.tab, tree.tip)
nmds.case.SparseDOSSA = get.unifrac(case.SparseDOSSA.phylo)

# (5) NorTA (assuming pre-computed)
case.NorTA.tab = otu_table(t(case.compo.NorTA.spec), taxa_are_rows = T)
case.NorTA.phylo = phyloseq(case.NorTA.tab, tree.tip)
nmds.case.NorTA = get.unifrac(case.NorTA.phylo)


#### 2. Process CONTROL group data (same order) ####
# (1) MBDDPM
ctrl.MBDDPM.filter = ctrl.MBDDPM[, colnames(ctrl.MBDDPM) %in% taxa.ctrl.NorTA]
ctrl.MBDDPM.tab = otu_table(t(ctrl.MBDDPM.filter), taxa_are_rows = T)
ctrl.MBDDPM.phylo = phyloseq(ctrl.MBDDPM.tab, tree.tip)
nmds.ctrl.MBDDPM = get.unifrac(ctrl.MBDDPM.phylo)

#  MBDDPM_NONE
ctrl.MBDDPM_NONE.filter = ctrl.MBDDPM_NONE[, colnames(ctrl.MBDDPM_NONE) %in% taxa.ctrl.NorTA]
ctrl.MBDDPM_NONE.tab = otu_table(t(ctrl.MBDDPM_NONE.filter), taxa_are_rows = T)
ctrl.MBDDPM_NONE.phylo = phyloseq(ctrl.MBDDPM_NONE.tab, tree.tip)
nmds.ctrl.MBDDPM_NONE = get.unifrac(ctrl.MBDDPM_NONE.phylo)

# (2) WGAN
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]
ctrl.WGAN.tab = otu_table(t(ctrl.WGAN.filter), taxa_are_rows = T)
ctrl.WGAN.phylo = phyloseq(ctrl.WGAN.tab, tree.tip)
nmds.ctrl.WGAN = get.unifrac(ctrl.WGAN.phylo)

# (3) MIDASim
ctrl.MIDASim.filter = ctrl.MIDASim[, colnames(ctrl.MIDASim) %in% taxa.ctrl.NorTA]
ctrl.MIDASim.tab = otu_table(t(ctrl.MIDASim.filter), taxa_are_rows = T)
ctrl.MIDASim.phylo = phyloseq(ctrl.MIDASim.tab, tree.tip)
nmds.ctrl.MIDASim = get.unifrac(ctrl.MIDASim.phylo)

# (4) SparseDOSSA2
ctrl.SparseDOSSA.filter = ctrl.compo.SparseDOSSA2[, colnames(ctrl.compo.SparseDOSSA2) %in% taxa.ctrl.NorTA]
ctrl.SparseDOSSA.tab = otu_table(t(ctrl.SparseDOSSA.filter), taxa_are_rows = T)
ctrl.SparseDOSSA.phylo = phyloseq(ctrl.SparseDOSSA.tab, tree.tip)
nmds.ctrl.SparseDOSSA = get.unifrac(ctrl.SparseDOSSA.phylo)

# (5) NorTA
ctrl.NorTA.tab = otu_table(t(ctrl.compo.NorTA.spec), taxa_are_rows = T)
ctrl.NorTA.phylo = phyloseq(ctrl.NorTA.tab, tree.tip)
nmds.ctrl.NorTA = get.unifrac(ctrl.NorTA.phylo)



#### 3. Save all results####
save.nam = paste0("./unifrac_nmds.Rdata")
save(
  # CONTROL group
  nmds.ctrl.real,
  nmds.ctrl.MBDDPM,
  nmds.ctrl.MBDDPM_NONE,
  nmds.ctrl.WGAN,
  nmds.ctrl.MIDASim,
  nmds.ctrl.SparseDOSSA,
  nmds.ctrl.NorTA,
  
  # CASE group
  nmds.case.real,
  nmds.case.MBDDPM,
  nmds.case.MBDDPM_NONE,
  nmds.case.WGAN,
  nmds.case.MIDASim,
  nmds.case.SparseDOSSA,
  nmds.case.NorTA,
  
  file = save.nam
)