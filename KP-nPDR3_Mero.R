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

setwd("/path/")

rep = c(1,2,3)
pos = c(0.251, 0, 0.45)
trf = c(0.051, 0, 0.03)
mer = c(0.148, 0.347, 0.926)
m_t = c(0.01, 0, 0.051)

Mero = data.frame(rep, pos, trf, mer, m_t)

Mero_long <- Mero  %>% gather("source", "value", 2:5)

Mero_long$source = as.factor(Mero_long$source)


Stat = ddply(Mero_long, c("source"), summarise, n = nrow(Mero),
             median=median(value),
             X25 = quantile(value, probs = 0.25), X75 = quantile(value, probs = 0.75),
             mean=mean(value),
             sd = sd(value),
             min = min(value), max = max(value))


names = c("1/3 MIC Meropenem + transferrin 750 µL/mL",
          "1/3 MIC Meropenem", "no drug", "transferrin 750 µL/mL")
Stat$source = names
write.table(Stat, file = "KP-nPDR#3_Mero_Statistics.txt", sep="\t", row.names=FALSE, quote=FALSE)
pdf("KP-nPDR#3_Mero_Statistics.pdf", height=22, width=12)
grid.table(Stat[1:9])
dev.off()

# Plot results
M = ggplot(data=Stat, aes(x=source, y=median)) +
  geom_bar(stat="identity", fill="deeppink3", width=0.5)+
  scale_x_discrete(name = NULL, labels = function(x) str_wrap(x, width = 10),
                   limits=c("no drug", "transferrin 750 µL/mL", "1/3 MIC Meropenem", 
                            "1/3 MIC Meropenem + transferrin 750 µL/mL")) +
  scale_y_continuous(name="# resistant mutants per 10^8 CFUs", limits = c(0, 0.7)) +
  ggtitle('Meropenem mutants selection in 24 hours. KP-nPDR#3') +
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
M

# Calculate stastic significance
wilcox.test(Mero$mer, Mero$m_t, paired=FALSE)
fit <- aov(value ~ source, data=Mero_long)
summary(fit) # display Type I ANOVA table
TukeyHSD(fit) # where fit comes from aov()
pvalues = TukeyHSD(fit)


# Save plot
ggsave("KP-nPDR#3_Mero_24h.pdf", width=10, height=7, dpi=200)
