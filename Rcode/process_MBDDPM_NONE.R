setwd("./data")
# summary data simulated by MB-DDPM_NONE #
cutoff = 1e-4

# load the MB-DDPM_NONE samples 
# 1) case group #
DDPM_NONE.case = read.table("epoch_150000_case_no-attention_04212100.csv", 
                       sep = ',', header = T)

# trim the small values #
DDPM_NONE.case[DDPM_NONE.case < cutoff] = 0


# 2) control group #
DDPM_NONE.ctrl = read.table("epoch_150000_ctrl_no-attention_04220920.csv",
                       sep = ',', header = T)

# # trim the small values #
DDPM_NONE.ctrl[DDPM_NONE.ctrl < cutoff] = 0


# rename all the taxa to the lowest rank #
name.case.ref = colnames(DDPM_NONE.case)
name.case.ref = sub('.*\\.', '',name.case.ref)
colnames(DDPM_NONE.case) = name.case.ref

name.ctrl.ref = colnames(DDPM_NONE.ctrl)
name.ctrl.ref = sub('.*\\.', '',name.ctrl.ref)
colnames(DDPM_NONE.ctrl) = name.ctrl.ref

rowSums(DDPM_NONE.case)
rowSums(DDPM_NONE.ctrl)

# convert to compositional data #
case.DDPM_NONE = DDPM_NONE.case; 
ctrl.DDPM_NONE = DDPM_NONE.ctrl

# save the results #
save(case.DDPM_NONE , ctrl.DDPM_NONE, file = "Data_by_MBDDPM_NONE_v1.0.Rdata")
name.case.DDPM_NONE = name.case.ref; name.ctrl.DDPM_NONE = name.ctrl.ref
save(name.case.DDPM_NONE, name.ctrl.DDPM_NONE, file = "taxa_name_MBDDPM_NONE_v1.0.Rdata")

