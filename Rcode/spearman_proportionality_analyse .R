setwd("D:/common/desktop/study/MB-DDPM/Rcode")
library(corrplot)
setwd("D:/common/desktop/study/MB-DDPM/data")
# load the names by NorTA (this is a subset)
load("taxa_name_NorTA.Rdata", verbose = T)

#### case group ####
# (a) load the real data #
load("Data_NielsenHB_real.Rdata", verbose = T)
load("Data_by_MBDDPM_v1.0.Rdata", verbose = T)
load("Data_by_WGAN_v1.0.Rdata", verbose = T)
load("Data_by_MIDASim_v1.0.Rdata", verbose = T)
load("Data_by_SparseDOSSA2_v1.0.Rdata", verbose = T)
load("Data_by_NorTA.Rdata", verbose = T)

# make sure the abundance matrices have the same spcies #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]
case.MBDDPM.filter = case.MBDDPM[, colnames(case.MBDDPM) %in% taxa.case.NorTA]
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]
case.MIDASim.filter = case.MIDASim[, colnames(case.MIDASim) %in% taxa.case.NorTA]
case.SparseDOSSA.filter = case.compo.SparseDOSSA2[, colnames(case.compo.SparseDOSSA2) %in% taxa.case.NorTA]
case.norta.filter = case.compo.NorTA.spec
# case.metap.filter = case.metap.full[, colnames(case.metap.full)%in% taxa.case.NorTA]
# set.seed(2147213)
# case.metap.filter = case.metap.filter[sample(1:nrow(case.metap.filter),size = 1000), ]

# extract the taxa to be included in the plot #
zero.real = apply(case.real.filter, 2, function(x) {sum(x == 0)}) == nrow(case.real.filter)
zero.DDPM = apply(case.MBDDPM.filter, 2, function(x) {sum(x == 0)}) == nrow(case.MBDDPM.filter)
zero.WGAN = apply(case.WGAN.filter, 2, function(x) {sum(x == 0)}) == nrow(case.WGAN.filter)
zero.MIDASim = apply(case.MIDASim.filter, 2, function(x) {sum(x == 0)}) == nrow(case.MIDASim.filter)
zero.SparseDOSSA = apply(case.SparseDOSSA.filter, 2, function(x) {sum(x == 0)}) == nrow(case.SparseDOSSA.filter)
zero.norta = apply(case.norta.filter, 2, function(x) {sum(x == 0)}) == nrow(case.norta.filter)
# zero.metap = apply(case.metap.filter, 2, function(x){sum(x == 0)}) == nrow(case.metap.filter)

threshold = 0.1
zero.count.by.taxa = apply(case.real.filter, 2, function(x) {sum(x == 0)})
extr.idx = zero.count.by.taxa < (nrow(case.real.filter) * threshold)

# get the name of taxa to be compared #
zero.rm.idx = zero.real | zero.DDPM | zero.WGAN | zero.MIDASim | 
  zero.SparseDOSSA | zero.norta
keep.idx = extr.idx & !zero.rm.idx
keep.nam = colnames(case.real.filter)[keep.idx]

# spearman correlation #
vis.spearman(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.WGAN.filter,
  Xmat4 = case.MIDASim.filter,
  Xmat5 = case.SparseDOSSA.filter,
  Xmat6 = case.norta.filter,
  label = "(a)",
  taxa.name = keep.nam
)

scatter.spearman(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.WGAN.filter,
  Xmat4 = case.MIDASim.filter,
  Xmat5 = case.SparseDOSSA.filter,
  Xmat6 = case.norta.filter,
  label = "(a)",
  taxa.name = keep.nam
)

# proportionality #
vis.proportionality(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.WGAN.filter,
  Xmat4 = case.MIDASim.filter,
  Xmat5 = case.SparseDOSSA.filter,
  Xmat6 = case.norta.filter,
  label = "(a)",
  taxa.name = keep.nam
)

scatter.proportionality(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.WGAN.filter,
  Xmat4 = case.MIDASim.filter,
  Xmat5 = case.SparseDOSSA.filter,
  Xmat6 = case.norta.filter,
  label = "(a)",
  taxa.name = keep.nam
)


#### 2. control group ####

# make sure the abundance matrices have the same spcies #
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
ctrl.MBDDPM.filter = ctrl.MBDDPM[, colnames(ctrl.MBDDPM) %in% taxa.ctrl.NorTA]
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]
ctrl.MIDASim.filter = ctrl.MIDASim[, colnames(ctrl.MIDASim) %in% taxa.ctrl.NorTA]
ctrl.SparseDOSSA.filter = ctrl.compo.SparseDOSSA2[, colnames(ctrl.compo.SparseDOSSA2) %in% taxa.ctrl.NorTA]
ctrl.norta.filter = ctrl.compo.NorTA.spec

# extract the taxa to be included in the plot #
# 1. remove the taxa with all zeros across all the samples #
zero.real = apply(ctrl.real.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.real.filter)
zero.DDPM = apply(ctrl.MBDDPM.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.MBDDPM.filter)
zero.WGAN = apply(ctrl.WGAN.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.WGAN.filter)
zero.MIDASim = apply(ctrl.MIDASim.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.MIDASim.filter)
zero.SparseDOSSA = apply(ctrl.SparseDOSSA.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.SparseDOSSA.filter)
zero.norta = apply(ctrl.norta.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.norta.filter)
# zero.metap = apply(ctrl.metap.filter, 2, function(x){sum(x == 0)}) == nrow(ctrl.metap.filter)

# 2. top abundant taxa in the real data #
threshold = 0.1
zero.count.by.taxa = apply(ctrl.real.filter, 2, function(x) {sum(x == 0)})
extr.idx = zero.count.by.taxa < (nrow(ctrl.real.filter) * threshold)

# 3. get the name of taxa to be compared #
zero.rm.idx = zero.real | zero.DDPM | zero.WGAN | zero.MIDASim | 
  zero.SparseDOSSA | zero.norta
keep.idx = extr.idx & !zero.rm.idx
keep.nam = colnames(ctrl.real.filter)[keep.idx]

# spearman correlation #
vis.spearman(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.WGAN.filter,
  Xmat4 = ctrl.MIDASim.filter,
  Xmat5 = ctrl.SparseDOSSA.filter,
  Xmat6 = ctrl.norta.filter,
  label = "(b)",
  taxa.name = keep.nam
)

scatter.spearman(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.WGAN.filter,
  Xmat4 = ctrl.MIDASim.filter,
  Xmat5 = ctrl.SparseDOSSA.filter,
  Xmat6 = ctrl.norta.filter,
  label = "(b)",
  taxa.name = keep.nam
)

# proportionality #
vis.proportionality(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.WGAN.filter,
  Xmat4 = ctrl.MIDASim.filter,
  Xmat5 = ctrl.SparseDOSSA.filter,
  Xmat6 = ctrl.norta.filter,
  label = "(b)",
  taxa.name = keep.nam
)

scatter.proportionality(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.WGAN.filter,
  Xmat4 = ctrl.MIDASim.filter,
  Xmat5 = ctrl.SparseDOSSA.filter,
  Xmat6 = ctrl.norta.filter,
  label = "(b)",
  taxa.name = keep.nam
)