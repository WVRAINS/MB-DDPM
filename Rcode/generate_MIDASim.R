# library(devtools)
# devtools::install_github("mengyu-he/MIDASim", build_vignettes = T, force = T)
# BiocManager::install("HMP2Data")
library(HMP2Data)
library(MIDASim)
library(scam)
library(future.apply)

setwd("./data")

load("NielsenHB_original_data.Rdata", verbose = TRUE)

species <- t(NielsenHB_count)

count.case <- species[patient.info[["disease"]] == "IBD", ]
count.ctrl <- species[patient.info[["disease"]] == "healthy", ]

case.setup <- MIDASim.setup(count.case, mode = 'parametric')
ctrl.setup <- MIDASim.setup(count.ctrl, mode = 'parametric')

n_case_samples <- 1000
n_ctrl_samples <- 1000

modified.case <-  MIDASim.modify(case.setup, 
                                 lib.size = sample(5000:20000, n_case_samples, replace = TRUE),
                                 mean.rel.abund = NULL,
                                 gengamma.mu = NULL,
                                 sample.1.prop = NULL,
                                 taxa.1.prop = NULL)
modified.ctrl <- MIDASim.modify(ctrl.setup,
                                lib.size = sample(5000:20000, n_ctrl_samples, replace = TRUE),
                                mean.rel.abund = NULL,
                                gengamma.mu = NULL,
                                sample.1.prop = NULL,
                                taxa.1.prop = NULL)

sim.case <- MIDASim(modified.case)
sim.ctrl <- MIDASim(modified.ctrl)

case.MIDASim <- sim.case$sim_rel
ctrl.MIDASim <- sim.ctrl$sim_rel

if (any(grepl("s__", colnames(case.MIDASim)))) {
  case_spec <- case.MIDASim[, grep("s__", colnames(case.MIDASim)), drop = FALSE]
  case.MIDASim <- t(apply(case_spec, 1, function(x) x / sum(x)))
} else {
  case.MIDASim <- t(apply(case.MIDASim, 1, function(x) x / sum(x)))
}

if (any(grepl("s__", colnames(ctrl.MIDASim)))) {
  ctrl_spec <- ctrl.MIDASim[, grep("s__", colnames(ctrl.MIDASim)), drop = FALSE]
  ctrl.MIDASim <- t(apply(ctrl_spec, 1, function(x) x / sum(x)))
} else {
  ctrl.MIDASim <- t(apply(ctrl.MIDASim, 1, function(x) x / sum(x)))
}

save(case.MIDASim, ctrl.MIDASim, file = "Data_by_MIDASim_v1.0.Rdata")