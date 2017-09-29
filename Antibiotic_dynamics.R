library(plyr)
library(readxl)


setwd("~/path")
df = read_excel("20d_R.xlsx", sheet = 1, col_names = TRUE, col_types = NULL, skip = 0)

pas <- df  %>% gather("source", "value", 3:6)

pas$strain = as.factor(pas$strain)
pas$source = as.factor(pas$source)

# subset data by antibiotic
list <- list()
strains = c("VA-AB21", "AB074", "PDR3", "PDR4")
for (i in strains){
  f = subset(pas, strain == i, select=day:value)
  list[[i]] <- f
}

plot_list = list()
labels = c("Ciprofloxacin", "Ciprofloxacin + transferrin", "Meropenem", "Meropenem + transferrin")

for (i in strains){
    p = ggplot(list[[i]], aes(colour = list[[i]]$source)) +
        geom_step(mapping=aes(x=list[[i]]$day, y=list[[i]]$value), size = 1.5) +
        geom_point(mapping=aes(x=list[[i]]$day, y=list[[i]]$value), size = 4) +
        ggtitle(i) +
        scale_x_continuous(name="Days", limits=c(0.5, 20.5), breaks = seq(0,20, 1)) +
        scale_y_log10(name="Antibiotic concentration",
                    breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)), limits = c(0.005, 50)) +
        annotation_logticks(sides = "l") +
        theme(panel.grid.major = element_line(linetype = "dotted", colour="grey26"), 
            panel.grid.major.x = element_blank(),
            axis.title.x = element_text(colour="grey26",size=20),
            axis.title.y = element_text(colour="grey26",size=20),
            plot.title = element_text(colour="grey26", size = 22, face = "bold.italic", vjust=0.5),
            legend.key.size = unit(1.4, "cm"),
            legend.key = element_rect(fill = "white", colour = "white"),
            legend.text = element_text(colour="grey26", size = 18),
            legend.title = element_text(colour="grey26", size = 22),
            axis.text = element_text(colour="grey26", size = 18),
            axis.line.x = element_line(colour = "grey26", size = 0.8),
            axis.line.y = element_line(colour = "grey26", size = 0.8),
            panel.background = element_blank()) +
            scale_colour_manual(name = "Treatment groups", labels = labels,
                                values = c("deeppink1", "deeppink4", "lightskyblue1", "skyblue3")) +
            scale_shape_manual(name = "Treatment groups", labels = labels, values = c(17, 19, 17, 19)) +
            geom_hline(yintercept=16, linetype = 2, color = "chartreuse3", size = 1.3) +
            annotate("text", 5, 20, label = "Susceptibility cutoff", size = 5)
        p
    plot_list[[i]] = p
}

pdf("plots_dynamics.pdf", width=13, height=7)
for (i in 1:4) {
  print(plot_list[[i]])
}
dev.off()


# Calculate non-parametric statistics
cipro_stat = list()
for (i in strains){
  a = subset(pas, strain == i, select=strain:value)
  a = subset(a, source =="cipro" | source =="cipro_trf", select=strain:value)
  cipro_stat[[i]] = wilcox.test(a$value ~ a$source, paired = TRUE)
}

save(cipro_stat, file = "cipro_stat.RData")

mero_stat = list()
for (i in strains){
  a = subset(pas, strain == i, select=strain:value)
  a = subset(a, source =="mero" | source =="mero_trf", select=strain:value)
  mero_stat[[i]] = wilcox.test(a$value ~ a$source, paired = TRUE)
}
save(mero_stat, file="mero_stat.RData")



