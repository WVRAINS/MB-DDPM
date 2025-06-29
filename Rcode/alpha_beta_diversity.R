setwd("D:/common/desktop/study/MB-DDPM/Rcode")
source("utility.R")
library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())
library(vegan)
library(energy)
library(entropy)
library(ggplot2)
library(cowplot)

library(scales) 
library(grDevices) 

setwd("D:/common/desktop/study/Diffusion-12.04-A/data")
load("taxa_name_NorTA.Rdata")

#### 1. case group ####

#### sample sparsity ####
load("Data_NielsenHB_real.Rdata", verbose = T)
case.real.filter = case.real.spec[, colnames(case.real.spec) %in% taxa.case.NorTA]

load("Data_by_MBDDPM_v1.0.Rdata")
case.MBDDPM.filter = case.MBDDPM[, colnames(case.MBDDPM) %in% taxa.case.NorTA]

load("Data_by_WGAN_v1.0.Rdata")
case.WGAN.filter = case.WGAN[, colnames(case.WGAN) %in% taxa.case.NorTA]

load("Data_by_MIDASim_v1.1.Rdata")
case.MIDASim.filter = case.MIDASim[, colnames(case.MIDASim) %in% taxa.case.NorTA]

load("Data_by_SparseDOSSA2_v1.0.Rdata")
case.SparseDOSSA.filter = case.compo.SparseDOSSA2[, colnames(case.compo.SparseDOSSA2) %in% taxa.case.NorTA]

load("Data_by_NorTA.Rdata")

shannon.case.real <- vegan::diversity(case.real.filter, index = "shannon")
shannon.case.MBDDPM <- vegan::diversity(case.MBDDPM.filter, index = "shannon")
shannon.case.WGAN <- vegan::diversity(case.WGAN.filter, index = "shannon")
shannon.case.MIDASim <- vegan::diversity(case.MIDASim.filter, index = "shannon")
shannon.case.SparseDOSSA <- vegan::diversity(case.SparseDOSSA.filter, index = "shannon")
shannon.case.NorTA <- vegan::diversity(case.compo.NorTA.spec, index = "shannon")


case.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"),
                                          c(length(shannon.case.real), length(shannon.case.MBDDPM),
                                            length(shannon.case.WGAN), length(shannon.case.MIDASim),
                                            length(shannon.case.SparseDOSSA), 1000)),
                               Diversity = c(shannon.case.real,
                                             shannon.case.MBDDPM,
                                             shannon.case.WGAN,
                                             shannon.case.MIDASim,
                                             shannon.case.SparseDOSSA,
                                             shannon.case.NorTA), stringsAsFactors = F)
case.diversity.df$Type <- factor(case.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"))
plot_colors <- c(
  "MB-DDPM" = "#FF0000",     
  "MB-GAN" = "#1E90FF",      
  "MIDASim" = "#2ca02c",     
  "SparseDOSSA" = "#BA55D3", 
  "NorTA" = "#FFD700",       
  "Real" = "#000000"         
)

case.diversity <- ggplot(case.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 7.5) + 
  ylab("Shannon Index") + scale_color_manual(values = plot_colors) +  
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),  
      c("Real", "MB-GAN"),
      c("Real", "MIDASim"),
      c("Real", "SparseDOSSA"),
      c("Real", "NorTA")
    ),label = "p.format"
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
simpsons.case.WGAN <- vegan::diversity(case.WGAN.filter, index = "simpson")
simpsons.case.MIDASim <- vegan::diversity(case.MIDASim.filter, index = "simpson")
simpsons.case.SparseDOSSA <- vegan::diversity(case.SparseDOSSA.filter, index = "simpson")
simpsons.case.NorTA <- vegan::diversity(case.compo.NorTA.spec, index = "simpson")

case.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"),
                                          c(length(simpsons.case.real), length(simpsons.case.MBDDPM),
                                            length(simpsons.case.WGAN), length(simpsons.case.MIDASim),
                                            length(simpsons.case.SparseDOSSA), 1000)),
                               Diversity = c(simpsons.case.real,
                                             simpsons.case.MBDDPM,
                                             simpsons.case.WGAN,
                                             simpsons.case.MIDASim,
                                             simpsons.case.SparseDOSSA,
                                             simpsons.case.NorTA), stringsAsFactors = F)
case.diversity.df$Type <- factor(case.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"))

case.diversity <- ggplot(case.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 1.6) + 
  ylab("Simpsons Index") + scale_color_manual(values = plot_colors) +  
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),  
      c("Real", "MB-GAN"),
      c("Real", "MIDASim"),
      c("Real", "SparseDOSSA"),
      c("Real", "NorTA")
    ),label = "p.format"
  )
case.diversity <- case.diversity + 
  theme(plot.margin = margin(b = 25))  
case.diversity <- ggdraw(case.diversity) + 
  draw_label(
    "(c)", 
    x = 0.02, y = 0.96,
    hjust = 0.5, vjust = 0, 
    size = 18,  
    fontface = "bold"  
  )
plot(case.diversity)


#### 2. contrl group ####
ctrl.real.filter = ctrl.real.spec[, colnames(ctrl.real.spec) %in% taxa.ctrl.NorTA]
ctrl.MBDDPM.filter = ctrl.MBDDPM[, colnames(ctrl.MBDDPM) %in% taxa.ctrl.NorTA]
ctrl.WGAN.filter = ctrl.WGAN[, colnames(ctrl.WGAN) %in% taxa.ctrl.NorTA]
ctrl.MIDASim.filter = ctrl.MIDASim[, colnames(ctrl.MIDASim) %in% taxa.ctrl.NorTA]
ctrl.SparseDOSSA.filter = ctrl.compo.SparseDOSSA2[, colnames(ctrl.compo.SparseDOSSA2) %in% taxa.ctrl.NorTA]

shannon.ctrl.real <- vegan::diversity(ctrl.real.filter, index = "shannon")
shannon.ctrl.MBDDPM <- vegan::diversity(ctrl.MBDDPM.filter, index = "shannon")
shannon.ctrl.WGAN <- vegan::diversity(ctrl.WGAN.filter, index = "shannon")
shannon.ctrl.MIDASim <- vegan::diversity(ctrl.MIDASim.filter, index = "shannon")
shannon.ctrl.SparseDOSSA <- vegan::diversity(ctrl.SparseDOSSA.filter, index = "shannon")
shannon.ctrl.NorTA <- vegan::diversity(ctrl.compo.NorTA.spec, index = "shannon")


ctrl.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"), 
                                          c(length(shannon.ctrl.real), length(shannon.ctrl.MBDDPM),
                                            length(shannon.ctrl.WGAN), length(shannon.ctrl.MIDASim),
                                            length(shannon.ctrl.SparseDOSSA), 1000)),
                               Diversity = c(shannon.ctrl.real,
                                             shannon.ctrl.MBDDPM,
                                             shannon.ctrl.WGAN,
                                             shannon.ctrl.MIDASim,
                                             shannon.ctrl.SparseDOSSA,
                                             shannon.ctrl.NorTA), stringsAsFactors = F)

ctrl.diversity.df$Type <- factor(ctrl.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"))

ctrl.diversity <- ggplot(ctrl.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 6.8) + 
  ylab("Shannon Index") + scale_color_manual(values = plot_colors) +  
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),  
      c("Real", "MB-GAN"),
      c("Real", "MIDASim"),
      c("Real", "SparseDOSSA"),
      c("Real", "NorTA")
    ),label = "p.format"
  )

ctrl.diversity <- ctrl.diversity + 
  theme(plot.margin = margin(b = 25))  
ctrl.diversity <- ggdraw(ctrl.diversity) + 
  draw_label(
    "(b)", 
    x = 0.02, y = 0.96,
    hjust = 0.5, vjust = 0, 
    size = 18,  
    fontface = "bold"  
  )
plot(ctrl.diversity)


simpsons.ctrl.real <- vegan::diversity(ctrl.real.filter, index = "simpson")
simpsons.ctrl.MBDDPM <- vegan::diversity(ctrl.MBDDPM.filter, index = "simpson")
simpsons.ctrl.WGAN <- vegan::diversity(ctrl.WGAN.filter, index = "simpson")
simpsons.ctrl.MIDASim <- vegan::diversity(ctrl.MIDASim.filter, index = "simpson")
simpsons.ctrl.SparseDOSSA <- vegan::diversity(ctrl.SparseDOSSA.filter, index = "simpson")
simpsons.ctrl.NorTA <- vegan::diversity(ctrl.compo.NorTA.spec, index = "simpson")

ctrl.diversity.df = data.frame(Type = rep(c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"), 
                                          c(length(simpsons.ctrl.real), length(simpsons.ctrl.MBDDPM),
                                            length(simpsons.ctrl.WGAN), length(simpsons.ctrl.MIDASim),
                                            length(simpsons.ctrl.SparseDOSSA), 1000)),
                               Diversity = c(simpsons.ctrl.real,
                                             simpsons.ctrl.MBDDPM,
                                             simpsons.ctrl.WGAN,
                                             simpsons.ctrl.MIDASim,
                                             simpsons.ctrl.SparseDOSSA,
                                             simpsons.ctrl.NorTA), stringsAsFactors = F)
ctrl.diversity.df$Type <- factor(ctrl.diversity.df$Type, 
                                 levels = c("Real", "MB-DDPM", "MB-GAN", "MIDASim", "SparseDOSSA", "NorTA"))

ctrl.diversity <- ggplot(ctrl.diversity.df, aes(x = Type, y = Diversity, color = Type)) +
  geom_boxplot(lwd = 0.8,fatten = 1.5,outlier.size = 1) + xlab(" ") +  ylim(0, 1.6) + 
  ylab("Simpsons Index") + scale_color_manual(values = plot_colors) +  
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  stat_compare_means(
    comparisons = list(
      c("Real", "MB-DDPM"),  
      c("Real", "MB-GAN"),
      c("Real", "MIDASim"),
      c("Real", "SparseDOSSA"),
      c("Real", "NorTA")
    ),label = "p.format"
  )
ctrl.diversity <- ctrl.diversity + 
  theme(plot.margin = margin(b = 25))  
ctrl.diversity <- ggdraw(ctrl.diversity) + 
  draw_label(
    "(d)", 
    x = 0.02, y = 0.96,
    hjust = 0.5, vjust = 0, 
    size = 18,  
    fontface = "bold"  
  )
plot(ctrl.diversity)



###    nMDS   ###
load("unifrac_nmds.Rdata", verbose = TRUE)

plot_nmds_comparison <- function(real, synthetic, synth_name, 
                                 color, main_title, xrange, yrange) {
  plot(1, 1, type = "n", main = main_title,
       xlim = xrange, ylim = yrange,
       xlab = "nMDS1", ylab = "nMDS2", 
       cex.lab = 1.5, cex.main = 1.8, cex.axis = 1.2)
  
  points(synthetic[, 1], synthetic[, 2], 
         pch = 20, col = scales::alpha(color, 0.5), cex = 1.2)
  points(real[, 1], real[, 2], 
         pch = 20, col = plot_colors["Real"], cex = 1.2)
  
  legend_x <- max(xrange) * 0.4
  legend_y <- max(yrange) * 0.99
  
  legend(x = legend_x, y = legend_y,
         legend = c("Real data", synth_name),
         col = c(plot_colors["Real"], color), 
         pch = 20, 
         pt.cex = 1.8,
         cex = 1.5,
         bty = "n",
         y.intersp = 1.2,
         x.intersp = 0.7,
         xpd = TRUE)
}

par(oma = c(2, 2, 2, 4))
par(mar = c(4, 4, 3, 1))

### Case group plots ###
layout(matrix(1:5, nrow = 1), widths = rep(1, 5))

ord.real <- nmds.case.real[["points"]]
ord_list <- list(
  "MB-DDPM" = nmds.case.MBDDPM[["points"]],
  "MB-GAN" = nmds.case.WGAN[["points"]],
  "MIDASim" = nmds.case.MIDASim[["points"]],
  "SparseDOSSA" = nmds.case.SparseDOSSA[["points"]],
  "NorTA" = nmds.case.NorTA[["points"]]
)

all_points <- c(list(ord.real), ord_list)
xrange <- range(sapply(all_points, function(x) range(x[,1]))) + c(-0.01, 0.01)
yrange <- range(sapply(all_points, function(x) range(x[,2]))) + c(-0.01, 0.01)

for (i in seq_along(ord_list)) {
  plot_nmds_comparison(ord.real, ord_list[[i]], names(ord_list)[i], 
                       plot_colors[names(ord_list)[i]], 
                       paste("Real vs", names(ord_list)[i]), 
                       xrange, yrange)
}

mtext("(a)", side = 3, line = -2, outer = TRUE, at = 0.01, cex = 1.8, font = 2)

### Control group plots ###
par(oma = c(2, 2, 2, 4))
par(mar = c(4, 4, 3, 1))
layout(matrix(1:5, nrow = 1), widths = rep(1, 5))

ord.real <- nmds.ctrl.real[["points"]]
ord_list <- list(
  "MB-DDPM" = nmds.ctrl.MBDDPM[["points"]],
  "MB-GAN" = nmds.ctrl.WGAN[["points"]],
  "MIDASim" = nmds.ctrl.MIDASim[["points"]],
  "SparseDOSSA" = nmds.ctrl.SparseDOSSA[["points"]],
  "NorTA" = nmds.ctrl.NorTA[["points"]]
)

all_points <- c(list(ord.real), ord_list)
xrange <- range(sapply(all_points, function(x) range(x[,1]))) + c(-0.01, 0.01)
yrange <- range(sapply(all_points, function(x) range(x[,2]))) + c(-0.1, 0.1)

for (i in seq_along(ord_list)) {
  plot_nmds_comparison(ord.real, ord_list[[i]], names(ord_list)[i], 
                       plot_colors[names(ord_list)[i]], 
                       paste("Real vs", names(ord_list)[i]), 
                       xrange, yrange)
}

mtext("(b)", side = 3, line = -2, outer = TRUE, at = 0.01, cex = 1.8, font = 2)