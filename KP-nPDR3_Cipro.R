library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(gridExtra)
library (plyr)
library(ggthemes)
library(scales)
library(svglite)
library(stringr)
library(ggsignif)

setwd("/path/")

rep = c(1,2,3)
pos = c(1.416, 1.419, 1.243)
trf = c(0.599, 1.329, 0)
cip = c(1.692, 1.949, 1.85)
c_t = c(0.853, 0.607, 0)

Cipro = data.frame(rep, pos, trf, cip, c_t)

Cipro_long <- Cipro  %>% gather("source", "value", 2:5)

Cipro_long$source = as.factor(Cipro_long$source)


Stat = ddply(Cipro_long, c("source"), summarise, n = nrow(Cipro),
             median=median(value),
             X25 = quantile(value, probs = 0.25), X75 = quantile(value, probs = 0.75),
             mean=mean(value),
             sd = sd(value),
             min = min(value), max = max(value))


names = c("1/3 MIC ciprofloxacin + transferrin 750 µL/mL",
          "1/3 MIC ciprofloxacin", "no drug", "transferrin 750 µL/mL")
Stat$source = names

# Save table
write.table(Stat, file = "KP-nPDR#3_Cipro_Statistics.txt", sep="\t", row.names=FALSE, quote=FALSE)
pdf("KP-nPDR#3_Cipro_Statistics.pdf", height=22, width=12)
grid.table(Stat[1:9])
dev.off()


# Plot results
C = ggplot(data=Stat, aes(x=source, y=median)) +
  geom_bar(stat="identity", fill="skyblue3", width=0.5)+
  scale_x_discrete(name = NULL, labels = function(x) str_wrap(x, width = 10),
                   limits=c("no drug", "transferrin 750 µL/mL", "1/3 MIC ciprofloxacin", 
                            "1/3 MIC ciprofloxacin + transferrin 750 µL/mL")) +
  scale_y_continuous(name="# resistant mutants per 10^8 CFUs", limits = c(0, 2.2)) +
  ggtitle('Ciprofloxacin mutants selection in 24 hours. KP-nPDR#3') +
  geom_errorbar(aes(ymin=X25, ymax=X75), width=0.1, size=0.6, color="midnightblue") +
  theme(panel.grid.major = element_line(linetype = "dotted", colour="grey26"), 
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text(colour="grey26",size=20),
        axis.title.y = element_text(colour="grey26",size=20),
        plot.title = element_text(colour="grey26", size = 22, face = "bold.italic", vjust=0.5),
        axis.text = element_text(colour="grey26", size = 18),
        axis.line.x = element_line(colour = "grey26", size = 0.8),
        axis.line.y = element_line(colour = "grey26", size = 0.8),
        panel.background = element_blank())
C


# Calculate stastic significance
wilcox.test(Cipro$cip, Cipro$c_t, paired=FALSE)
fit <- aov(value ~ source, data=Cipro_long)
summary(fit) # display Type I ANOVA table
TukeyHSD(fit) # where fit comes from aov()
pvalues = TukeyHSD(fit)


# Add significance to the plot
C = C + geom_signif(stat ="identity", size=1, textsize=6,
                    aes(x=3,xend=4, y=2.1, yend=2.1, annotation= "p < 0.05")) +
    geom_signif(stat ="identity", size=1, textsize=6,
                    aes(x=2,xend=3, y=1.95, yend=1.95, annotation= "p < 0.05"))
C

# Save plot
ggsave("KP-nPDR#3_Cipro_24h.pdf", width=10, height=7, dpi=200)
