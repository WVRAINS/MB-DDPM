setwd("D:/common/desktop/study/MB-DDPM/Rcode")
source("utility.R")
library(corrplot)
library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())
library(phyloseq)
library(vegan)
library(energy)
library(entropy)
library(ggplot2)
library(cowplot)

setwd("D:/common/desktop/study/MB-DDPM/data")
load("taxa_name_NorTA.Rdata")

#### case group ####

#### sample sparsity ####
load("Data_NielsenHB_real.Rdata", verbose = T)
# restrict to the subgroup by case group #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]

load("Data_by_MBDDPM_v1.0.Rdata")
# restrict to the subgroup by case group #
case.MBDDPM.filter = case.MBDDPM[, colnames(case.MBDDPM) %in% taxa.case.NorTA]

load("Data_by_MBDDPM_NONE_v1.0.Rdata")
# restrict to the subgroup by case group #
case.MBDDPM_NONE.filter = case.MBDDPM_NONE[, colnames(case.MBDDPM_NONE) %in% taxa.case.NorTA]

#### 1.1 shannon index ####
shannon.case.real <- vegan::diversity(case.real.filter, index = "shannon")
shannon.case.MBDDPM <- vegan::diversity(case.MBDDPM.filter, index = "shannon")
shannon.case.MBDDPM_NONE <- vegan::diversity(case.MBDDPM_NONE.filter, index = "shannon")

case.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-DDPM_NONE"),
                                          c(length(shannon.case.real), length(shannon.case.MBDDPM),
                                            length(shannon.case.MBDDPM_NONE))),
                               Diversity = c(shannon.case.real,
                                             shannon.case.MBDDPM,
                                             shannon.case.MBDDPM_NONE), stringsAsFactors = F)
case.diversity.df$Type <- factor(case.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-DDPM_NONE"))
plot_colors <- c(
  "MB-DDPM" = "#FF0000", 
  "MB-DDPM_NONE" = "#1E90FF", 
  "Real" = "#000000" 
)

# 绘图
case.diversity <- ggplot(case.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 4.5) + 
  ylab("Shannon Index") + scale_color_manual(values = plot_colors) +  # 使用自定义颜色
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),
      c("Real", "MB-DDPM_NONE")
    ),label = "p.signif"
  )

case.diversity <- case.diversity + 
  theme(plot.margin = margin(b = 25)) 
case.diversity <- ggdraw(case.diversity) + 
  draw_label(
    "(a)", 
    x = 0.02, y = 0.98,
    size = 18, 
    fontface = "bold"
  )

plot(case.diversity)

simpsons.case.real <- vegan::diversity(case.real.filter, index = "simpson")
simpsons.case.MBDDPM <- vegan::diversity(case.MBDDPM.filter, index = "simpson")
simpsons.case.MBDDPM_NONE <- vegan::diversity(case.MBDDPM_NONE.filter, index = "simpson")

case.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-DDPM_NONE"),
                                          c(length(simpsons.case.real), length(simpsons.case.MBDDPM),
                                            length(simpsons.case.MBDDPM_NONE))),
                               Diversity = c(simpsons.case.real,
                                             simpsons.case.MBDDPM,
                                             simpsons.case.MBDDPM_NONE), stringsAsFactors = F)
case.diversity.df$Type <- factor(case.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-DDPM_NONE"))
plot_colors <- c(
  "MB-DDPM" = "#FF0000",
  "MB-DDPM_NONE" = "#1E90FF",
  "Real" = "#000000"
)

# 绘图
case.diversity <- ggplot(case.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0.5, 1.1) + 
  ylab("Simpson Index") + scale_color_manual(values = plot_colors) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),
      c("Real", "MB-DDPM_NONE")
    ),label = "p.signif"
  )

case.diversity <- case.diversity + 
  theme(plot.margin = margin(b = 25))
case.diversity <- ggdraw(case.diversity) + 
  draw_label(
    "(c)", 
    x = 0.02, y = 0.98,
    size = 18,
    fontface = "bold" 
  )
plot(case.diversity)


##### nMDS #####
# Load results
load("unifrac_nmds.Rdata", verbose = TRUE)

plot_colors <- c(
  "MB-DDPM" = "#FF0000",
  "MB-DDPM_NONE" = "#1E90FF",
  "Real" = "#000000"
)

plot_nmds_comparison <- function(real, synthetic, synth_name, 
                                 color, main_title, xrange, yrange) {
  plot(1, 1, type = "n", main = main_title,
       xlim = xrange, ylim = yrange,
       xlab = "nMDS1", ylab = "nMDS2", 
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2)
  
  points(synthetic[, 1], synthetic[, 2], 
         pch = 20, col = scales::alpha(color, 0.5), cex = 1.2)
  points(real[, 1], real[, 2], 
         pch = 20, col = plot_colors["Real"], cex = 1.2)
  
  legend(x = "topright",
         legend = c("Real data", synth_name),
         col = c(plot_colors["Real"], color), 
         pch = 20, 
         pt.cex = 1.8,
         cex = 1.3, 
         bty = "n", 
         y.intersp = 1.2,
         x.intersp = 1.2,
         inset = c(0, 0))
}

### Case ###
par(mar = c(4, 6, 2, 6))
par(mfcol = c(1, 2))

ord.real <- nmds.case.real[["points"]]
ord_list <- list(
  "MB-DDPM" = nmds.case.MBDDPM[["points"]],
  "MB-DDPM_NONE" = nmds.case.MBDDPM_NONE[["points"]]
)

all_points <- c(list(ord.real), ord_list)
xrange <- range(sapply(all_points, function(x) range(x[,1]))) + c(-0.01, 0.01)
yrange <- range(sapply(all_points, function(x) range(x[,2]))) + c(-0.01, 0.01)

for (i in seq_along(ord_list)) {
  plot_nmds_comparison(ord.real, ord_list[[i]], names(ord_list)[i], 
                       plot_colors[names(ord_list)[i]], paste("Real vs", names(ord_list)[i]), 
                       xrange, yrange)
}
mtext("(a)", side = 3, line = -2, outer = TRUE, at = 0.01, cex = 1.5, font = 2)

### ctrl ###
par(mar = c(4, 6, 2, 6))
par(mfcol = c(1, 2))

ord.real <- nmds.ctrl.real[["points"]]
ord_list <- list(
  "MB-DDPM" = nmds.ctrl.MBDDPM[["points"]],
  "MB-DDPM_NONE" = nmds.ctrl.MBDDPM_NONE[["points"]]
)

all_points <- c(list(ord.real), ord_list)
xrange <- range(sapply(all_points, function(x) range(x[,1]))) + c(-0.01, 0.01)
yrange <- range(sapply(all_points, function(x) range(x[,2]))) + c(-0.01, 0.01)

for (i in seq_along(ord_list)) {
  plot_nmds_comparison(ord.real, ord_list[[i]], names(ord_list)[i], 
                       plot_colors[names(ord_list)[i]], paste("Real vs", names(ord_list)[i]), 
                       xrange, yrange)
}
mtext("(b)", side = 3, line = -2, outer = TRUE, at = 0.01, cex = 1.5, font = 2)



#### ctrl group ####
# restrict to the subgroup by ctrl group #
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
# (b) MB-DDPM #
# restrict to the subgroup by ctrl group #
ctrl.MBDDPM.filter = ctrl.MBDDPM[, colnames(ctrl.MBDDPM) %in% taxa.ctrl.NorTA]
# restrict to the subgroup by ctrl group #
ctrl.MBDDPM_NONE.filter = ctrl.MBDDPM_NONE[, colnames(ctrl.MBDDPM_NONE) %in% taxa.ctrl.NorTA]

#### 1.1 shannon index ####
shannon.ctrl.real <- vegan::diversity(ctrl.real.filter, index = "shannon")
shannon.ctrl.MBDDPM <- vegan::diversity(ctrl.MBDDPM.filter, index = "shannon")
shannon.ctrl.MBDDPM_NONE <- vegan::diversity(ctrl.MBDDPM_NONE.filter, index = "shannon")

ctrl.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-DDPM_NONE"),
                                          c(length(shannon.ctrl.real), length(shannon.ctrl.MBDDPM),
                                            length(shannon.ctrl.MBDDPM_NONE))),
                               Diversity = c(shannon.ctrl.real,
                                             shannon.ctrl.MBDDPM,
                                             shannon.ctrl.MBDDPM_NONE), stringsAsFactors = F)
ctrl.diversity.df$Type <- factor(ctrl.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-DDPM_NONE"))
plot_colors <- c(
  "MB-DDPM" = "#FF0000",
  "MB-DDPM_NONE" = "#1E90FF",  
  "Real" = "#000000" 
)

# 绘图
ctrl.diversity <- ggplot(ctrl.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 4.5) + 
  ylab("Shannon Index") + scale_color_manual(values = plot_colors) +  # 使用自定义颜色
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),
      c("Real", "MB-DDPM_NONE")
    ),label = "p.signif"
  )

ctrl.diversity <- ctrl.diversity + 
  theme(plot.margin = margin(b = 25))
ctrl.diversity <- ggdraw(ctrl.diversity) + 
  draw_label(
    "(b)", 
    x = 0.02, y = 0.98,
    size = 18,
    fontface = "bold"
  )

plot(ctrl.diversity)

simpsons.ctrl.real <-  vegan::diversity(ctrl.real.filter, index = "simpson")
simpsons.ctrl.MBDDPM <-  vegan::diversity(ctrl.MBDDPM.filter, index = "simpson")
simpsons.ctrl.MBDDPM_NONE <-  vegan::diversity(ctrl.MBDDPM_NONE.filter, index = "simpson")

ctrl.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-DDPM_NONE"),
                                          c(length(simpsons.ctrl.real), length(simpsons.ctrl.MBDDPM),
                                            length(simpsons.ctrl.MBDDPM_NONE))),
                               Diversity = c(simpsons.ctrl.real,
                                             simpsons.ctrl.MBDDPM,
                                             simpsons.ctrl.MBDDPM_NONE), stringsAsFactors = F)
ctrl.diversity.df$Type <- factor(ctrl.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-DDPM_NONE"))
plot_colors <- c(
  "MB-DDPM" = "#FF0000",
  "MB-DDPM_NONE" = "#1E90FF", 
  "Real" = "#000000"
)

# 绘图
ctrl.diversity <- ggplot(ctrl.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 1.1) + 
  ylab("Simpson Index") + scale_color_manual(values = plot_colors) + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),
      c("Real", "MB-DDPM_NONE")
    ),label = "p.signif"
  )

ctrl.diversity <- ctrl.diversity + 
  theme(plot.margin = margin(b = 25)) 
ctrl.diversity <- ggdraw(ctrl.diversity) + 
  draw_label(
    "(d)", 
    x = 0.02, y = 0.98,
    size = 18, 
    fontface = "bold" 
  )

plot(ctrl.diversity)


#### case group ####
# make sure the abundance matrices have the same spcies #
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]
case.MBDDPM.filter = case.MBDDPM[, colnames(case.MBDDPM) %in% taxa.case.NorTA]
case.MBDDPM_NONE.filter = case.MBDDPM_NONE[, colnames(case.MBDDPM_NONE) %in% taxa.case.NorTA]

# extract the taxa to be included in the plot #
# remove the taxa with all zeros across all the samples
zero.real = apply(case.real.filter, 2, function(x) {sum(x == 0)}) == nrow(case.real.filter)
zero.DDPM = apply(case.MBDDPM.filter, 2, function(x) {sum(x == 0)}) == nrow(case.MBDDPM.filter)
zero.DDPM_NONE = apply(case.MBDDPM_NONE.filter, 2, function(x) {sum(x == 0)}) == nrow(case.MBDDPM_NONE.filter)

# top abundant taxa in the real data
threshold = 0.1
zero.count.by.taxa = apply(case.real.filter, 2, function(x) {sum(x == 0)})
extr.idx = zero.count.by.taxa < (nrow(case.real.filter) * threshold)

# get the name of taxa to be compared #
zero.rm.idx = zero.real | zero.DDPM | zero.DDPM_NONE
keep.idx = extr.idx & !zero.rm.idx
keep.nam = colnames(case.real.filter)[keep.idx]

# spearman correlation #
vis.spearman(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(a)",
  titles = c("Real Data", "MB-DDPM", "MB-DDPM_NONE")
)

scatter.spearman(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(a)",
  method.names = c("MB-DDPM", "MB-DDPM_NONE")
)

# proportionality #
vis.proportionality(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(a)",
  titles = c("Real Data", "MB-DDPM", "MB-DDPM_NONE")
)

scatter.proportionality(
  Xmat1 = case.real.filter,
  Xmat2 = case.MBDDPM.filter,
  Xmat3 = case.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(a)",
  method.names = c("MB-DDPM", "MB-DDPM_NONE")
)

# make sure the abundance matrices have the same spcies #
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
ctrl.MBDDPM.filter = ctrl.MBDDPM[, colnames(ctrl.MBDDPM) %in% taxa.ctrl.NorTA]
ctrl.MBDDPM_NONE.filter = ctrl.MBDDPM_NONE[, colnames(ctrl.MBDDPM_NONE) %in% taxa.ctrl.NorTA]

# extract the taxa to be included in the plot #
# remove the taxa with all zeros across all the samples
zero.real = apply(ctrl.real.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.real.filter)
zero.DDPM = apply(ctrl.MBDDPM.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.MBDDPM.filter)
zero.DDPM_NONE = apply(ctrl.MBDDPM_NONE.filter, 2, function(x) {sum(x == 0)}) == nrow(ctrl.MBDDPM_NONE.filter)

# top abundant taxa in the real data
threshold = 0.1
zero.count.by.taxa = apply(ctrl.real.filter, 2, function(x) {sum(x == 0)})
extr.idx = zero.count.by.taxa < (nrow(ctrl.real.filter) * threshold)

# get the name of taxa to be compared #
zero.rm.idx = zero.real | zero.DDPM | zero.DDPM_NONE
keep.idx = extr.idx & !zero.rm.idx
keep.nam = colnames(ctrl.real.filter)[keep.idx]

# spearman correlation #
vis.spearman(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(b)",
  titles = c("Real Data", "MB-DDPM", "MB-DDPM_NONE")
)

scatter.spearman(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(b)",
  method.names = c("MB-DDPM", "MB-DDPM_NONE")
)

# proportionality #
vis.proportionality(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(b)",
  titles = c("Real Data", "MB-DDPM", "MB-DDPM_NONE")
)

scatter.proportionality(
  Xmat1 = ctrl.real.filter,
  Xmat2 = ctrl.MBDDPM.filter,
  Xmat3 = ctrl.MBDDPM_NONE.filter,
  taxa.name = keep.nam,
  label = "(b)",
  method.names = c("MB-DDPM", "MB-DDPM_NONE")
)

