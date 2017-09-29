library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(gridExtra)
library (plyr)
library(ggthemes)
library(scales)
library(svglite)

# Export data of all experiments to the data frame

Plate_names = c("Time", "0_Cipro", "1/3MIC_trf+Cipro", 
                "1MIC_trf+Cipro", "3MIC_trf+Cipro",
                "0_Mero", "1/3MIC_trf+Mero",
                "1MIC_trf+Mero", "3MIC_trf+Mero",
                "1/3MIC_trf", "1MIC_trf", "3MIC_trf")

Exp = c(5:8)                                                      #CHANGE HERE

df <- list()
setwd("/Users/ksun/Dropbox/USC/Transferrin_experiments/Time-kill_Assay/UH118/1'3MIC")                         #CHANGE HERE
for (i in Exp){
  f = paste("Result_corrected_1'3MIC_exp", i, sep="_")
  f = paste(f,"txt", sep = ".")
  print(f)
  df[[i]] <- read.table(f, sep = "\t", header = FALSE, skip = 1)
  # print(head(df))
  colnames(df[[i]]) = Plate_names
}

# Calculate Reduction = (F-Fneg)/(Fpos-Fneg)
df1 = list()
for (i in Exp){
  test = df[[i]]
  vec = unname(as.numeric(test[1,]), force = FALSE)
  df1[[i]] = sweep(test, 2, vec, `/`)
}

# Gather data
for (i in Exp){
  df1[[i]] <- gather(df1[[i]], "source", "value", 2:length(df1[[i]])) %>% 
    separate("source", c("Dose_of_trf", "Group", "Replica"), sep= '_')
  df1[[i]]$Replica <- i
}

df2 <- df1[[1]]
for (i in Exp){
  df2 <- rbind(df2, df1[[i]])
}

# Replace Time column and remove posirive control
df2 = df2 %>% drop_na()
c = rep(c(1,2,4,6,8,24), 44)
df2$Time = c

write.table(df2, file = "All_data_long_UH118_1'3MIC.txt", sep="\t", row.names=FALSE, quote=FALSE)

# set up a group factors
df2$Dose_of_trf = factor(df2$Dose_of_trf, levels = unique(df2$Dose_of_trf))
df2$Group = factor(df2$Group, levels = unique(df2$Group))

# Calculate descriptive statistics
Stat = ddply(df2, c("Group", "Dose_of_trf", "Time"), summarise, n = length(unique(df2$Replica)), mean=mean(value),
             sd = sd(value), median=median(value),
             X25 = quantile(value, probs = 0.25), X75 = quantile(value, probs = 0.75),
             min = min(value), max = max(value))

# Write down summary tables
write.table(Stat, file = "Descriptive_statistics_UH118_1'3MIC.txt", sep="\t", row.names=FALSE, quote=FALSE)
pdf("Descriptive_statistics_UH118_1'3MIC.pdf", height=22, width=12)
grid.table(Stat[1:11])
dev.off()

# Add name column
Stat$name = paste(Stat$Dose_of_trf, Stat$Group, sep="_")

# subset data by antibiotic
Cipro <- subset(Stat, Group=="trf" | Group=="trf+Cipro" | Group=="Cipro",
                select=Group:name)
Mero <- subset(Stat, Group=="trf" | Group=="trf+Mero" | Group=="Mero",
               select=Group:name)

# Plot results; Cipro
labC = c("Cipro", "1/3 MIC trf", "1/3 MIC trf + 1/3 MIC Cipro", "1 MIC trf", 
         "1 MIC trf + 1/3 MIC Cipro", "3 MIC trf", "3 MIC trf + 1/3 MIC Cipro")

C = ggplot(Cipro, aes(x = Time, y = median, colour=name, shape = name)) +
  geom_line(size = 1.1) +
  geom_point(size = 4.5)
# change axises limits and add title
C = C + scale_x_continuous(name="Time (Hours)", limits=c(0.5, 24.3), breaks = seq(0,24, 4)) +
  scale_y_continuous(name="Resazurin reduction", limits = c(0, 1)) +
  ggtitle('A. baumannii UH118; 1/3x Ciprofloxacin MIC; Time-kill') +
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
        panel.background = element_blank())
C = C + scale_colour_manual(name = "Treatment Groups",
                            labels = labC,
                            values = c("black", "deeppink3", "skyblue3", "deeppink3", "skyblue3", 
                                       "deeppink3", "skyblue3")) +
  scale_shape_manual(name = "Treatment Groups",
                     labels = labC,
                     values = c(18, 17, 17, 19, 19, 15, 15))
C = C + geom_errorbar(aes(ymin=X25, ymax=X75), width=0.3, size=0.6, color="midnightblue")
C
# Save summary chart
ggsave("UH118_1'3MIC_Cipro_new.pdf", width=15, height=7, dpi=200)
ggsave("/Users/ksun/Dropbox/USC/Transferrin_experiments/Time-kill_Assay/Final_graphs/UH118_1'3MIC_Cipro_new.pdf", width=15, height=7, dpi=200)
ggsave(file = "/Users/ksun/Dropbox/USC/Transferrin_experiments/Time-kill_Assay/Final_graphs/svg/UH118_1'3MIC_Cipro_new.svg", C, width=20, height=9, dpi=400)

# Plot results; Mero
labM = c("Mero", "1/3 MIC trf", "1/3 MIC trf + 1/3 MIC Mero", "1 MIC trf", 
         "1 MIC trf + 1/3 MIC Mero", "3 MIC trf", "3 MIC trf + 1/3 MIC Mero")

M = ggplot(Mero, aes(x = Time, y = median, colour=name, shape = name)) +
  geom_line(size = 1.1) +
  geom_point(size = 4.5)
# change axises limits and add title
M = M + scale_x_continuous(name="Time (Hours)", limits=c(0.5, 24.3), breaks = seq(0,24, 4)) +
  scale_y_continuous(name="Resazurin reduction", limits = c(0, 1)) +
  ggtitle('A. baumannii UH118; 1/3x Meropenem MIC; Time Kill') +
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
        panel.background = element_blank())
M = M + scale_colour_manual(name = "Treatment Groups",
                            labels = labM,
                            values = c("black", "deeppink3", "skyblue3", "deeppink3", "skyblue3", 
                                       "deeppink3", "skyblue3")) +
  scale_shape_manual(name = "Treatment Groups",
                     labels = labM,
                     values = c(18, 17, 17, 19, 19, 15, 15))
M = M + geom_errorbar(aes(ymin=X25, ymax=X75), width=0.3, size=0.6, color="midnightblue")
M
# Save summary chart
ggsave("UH118_1'3MIC_Mero_new.pdf", width=15, height=7, dpi=200)
ggsave("/Users/ksun/Dropbox/USC/Transferrin_experiments/Time-kill_Assay/Final_graphs/UH118_1'3MIC_Mero_new.pdf", width=15, height=7, dpi=200)
ggsave(file = "/Users/ksun/Dropbox/USC/Transferrin_experiments/Time-kill_Assay/Final_graphs/svg/UH118_1'3MIC_Mero_new.svg", M, width=20, height=9, dpi=400)

