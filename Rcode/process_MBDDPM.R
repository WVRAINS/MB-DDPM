setwd("./data")
# summary data simulated by MB-DDPM #
cutoff = 1e-4

# load the MB-DDPM samples 
# 1) case group #
DDPM.case = read.table("epoch_150000_case_fold_04202100.csv", 
                      sep = ',', header = T)

# trim the small values #
DDPM.case[DDPM.case < cutoff] = 0


# 2) control group #
DDPM.ctrl = read.table("epoch_150000_ctrl_fold_04172340.csv",
                      sep = ',', header = T)

# # trim the small values #
DDPM.ctrl[DDPM.ctrl < cutoff] = 0


# rename all the taxa to the lowest rank #
name.case.ref = colnames(DDPM.case)
name.case.ref = sub('.*\\.', '',name.case.ref)
colnames(DDPM.case) = name.case.ref

name.ctrl.ref = colnames(DDPM.ctrl)
name.ctrl.ref = sub('.*\\.', '',name.ctrl.ref)
colnames(DDPM.ctrl) = name.ctrl.ref

rowSums(DDPM.case)
rowSums(DDPM.ctrl)

# convert to compositional data #
case.MBDDPM = DDPM.case; 
ctrl.MBDDPM = DDPM.ctrl

# save the results #
save(case.MBDDPM , ctrl.MBDDPM, file = "Data_by_MBDDPM_v1.0.Rdata")
name.case.MBDDPM = name.case.ref; name.ctrl.MBDDPM = name.ctrl.ref
save(name.case.MBDDPM, name.ctrl.MBDDPM, file = "taxa_name_DDPM_v1.0.Rdata")

