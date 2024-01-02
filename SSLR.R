source("./ggplot2_config.R")
library(ggalluvial)
library(reshape2)
library(export)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(dplyr)
library(plyr) 
library(phyloseq)
.libPaths("./R-4.2.0/library")
.libPaths("./win-library")
#Fig1

setwd("./DNA_Mock")

TPM <- read.csv("./DNA_Mock/FV_result/Mock_percentage_final.csv",sep = ',', header = T)

TPM$Method = factor(TPM$Method, levels=c('ExpectedA', 'ExpectedB','ExpectedC','SSLR10','SSLR100','xGen','MDA_0.5h','MDA_1.5h','Nextera'))
TPM$Phage = factor(TPM$Phage, levels=c('T4','T7','P1','Lambda','C2','P35','Phi29','PhiX174','M13mp18'))

TPM1 <- TPM %>% filter(Method %in% c('ExpectedA', 'ExpectedB','ExpectedC')) 
Fig2A_DNA <- TPM1 %>% filter(c(T4=="Yes")) %>%
  ggplot(aes(x= Method,y=real,fill=Phage))+
  geom_bar(stat="identity",position="stack")+
  labs(x="", y="Percentage of genomes (%)") +
  facet_grid(~Mock, scales= "free_x") +
  scale_fill_manual(values =my_palette_large[1:10])+
  theme_classic() +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10),
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1)) +
  ggtitle("DNA phage Mock without T4_expected")
graph2ppt(Fig2A_DNA, file="./DNA_Mock/FV_result/DNA phage Mock without T4_expected",paper= "A4", width = 5, height = 5, scaling = 75, vector.graphic = TRUE)


#Read in raw otu data 
otu_mat1 <- read.delim("raw_vOTU_table.txt",sep="\t",header=T)
rownames(otu_mat1) <- otu_mat1[,1]
otu_mat1 <- otu_mat1[,-1]
otu_mat1 <- as.matrix(otu_mat1)


#Read in taxonomy data
tax_mat <- read.table("Taxonomy_viral.csv",sep = ',',header = T)
tax_mat <- as.matrix(tax_mat)
rownames(tax_mat)<-tax_mat[,1]
tax_mat<-tax_mat[,-1]
TAX=tax_table(tax_mat)

#Read in metadata
sample_df <- read.delim("metadata.txt",sep="\t",header=T)
row.names(sample_df) <- sample_df[,1] 
#sample_df <- sample_df[,-1]

SAM <- sample_data(sample_df,errorIfNULL = T)
SAM1 <- SAM %>% filter(!Method %in% c('SSLR100')) 
OTU <- otu_table(otu_mat1, taxa_are_rows = TRUE)


ps <- phyloseq(OTU, TAX, SAM)


ps <- subset_samples(ps, Method!="SSLR100")

my_palette_large <- c("#00008B","#FFB90F","#8FBC8F","#9932CC","#CAFF70","#E64B35FF")
ps@sam_data$Method = factor(ps@sam_data$Method, levels=c("SSLR10","xGen","MDA_0.5h","MDA_1.5h","Nextera"))
#ps@sam_data$Method = factor(ps@sam_data$Method, levels=c("SSLR","xGen","MDA_0.5h","MDA_1.5h","Nextera"))
dist = phyloseq::distance(ps, method="bray", binary = FALSE)
ordination = ordinate(ps, method="PCoA", distance=dist)
p11 <- plot_ordination(ps, ordination, color="Method", shape = "Mock") +
  geom_point(size=3) +
  scale_fill_manual(values = my_palette_large) +
  scale_fill_manual(values = my_palette_large) +
  geom_vline(xintercept = 0,lty=2,col="grey50") +
  geom_hline(yintercept = 0,lty=2,col="grey50") +
  geom_text(label = ps@sam_data$T4, size = 5, hjust = 0, nudge_x = 0.0025) +
  #geom_label_repel(aes(ps, ordination, fill=factor(Method), label=SAM$Method)) + 
  theme_classic() +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10),
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1))
p11
graph2ppt(p11, file="FV_result/bray_pcoA1",paper= "A4", width = 7, height = 5, scaling = 50, vector.graphic = TRUE)
##################  dsDNA vs ssDNA ##################################################
#Lambda1 <- melt(Lambda)
Fig1_lambda_DNA <- ggplot(Lambda,aes(x= Method,y=Relative.abundance...., fill= Phage))+
  geom_bar(stat="identity",position="stack")+
  labs(x="", y="Abundance level (%)") +
  facet_wrap(~CA, scales= "free_x",nrow =1) +
  scale_fill_manual(values =my_palette_large[22:29])+
  geom_text(aes(label = Relative.abundance....), position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 14),#设定分页图的标题属性
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10),
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1)) +
  ggtitle("Expected phage abundance level")
graph2ppt(Fig1_lambda_DNA, file="Fig_dsDNA_ssDNA_remove_T4",paper= "A4",width = 15, height = 10,vector.graphic = TRUE)


Lambda <- read.csv("Fig2_dsDNA.csv",sep = ',', header = T)
Lambda1 <- melt(Lambda)
Lambda1$variable = factor(Lambda1$variable, levels=c('Expected', 'EA1_SSLR', 'EA2_SSLR',
                                               'EB1_SSLR','EB2_SSLR', 
                                                'EC1_SSLR','EC2_SSLR'))

Fig1_lambda_DNA <- ggplot(Lambda1,aes(x=variable,y=value, fill= Phage))+
  geom_bar(stat="identity",position="stack")+
  labs(x="", y="Abundance level (%)") +
  
  scale_fill_manual(values =my_palette_large[22:29])+
  theme_classic() +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 14),#设定分页图的标题属性
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10),
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1)) +
  ggtitle("Expected lambda phage abundance level")
graph2ppt(Fig1_lambda_DNA, file="Fig3D_dsDNA_SSLR",paper= "A4",width = 7.5, height = 10,vector.graphic = TRUE)

TPM1 <- ddply(TPM,~ TPM,percent = 1/sum(TPM)*100)
ggplot(data,aes(stage,percent,fill=gender))+
  geom_bar(stat="identity",position="stack")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))


x1<- MDA %>% filter(V1 == c("NC_003714.1"))
x2 <- sslr %>% filter(V1 == c("NC_003714.1"))
write.table(x2,"x2.txt",sep="\t",quote=F,row.names=T)

Phi6 <- read.table("x.txt", sep='\t', header=T)
Phi6 <- as.data.frame(Phi6)
Phi6 <- melt(Phi6,id.vars = c("Phi6_S","Position"),measure.vars = c("MDA","SSLR"))

Phi6 <- ggplot(Phi6,aes(x=Position,y=value,color=variable)) +
  geom_line() +
  scale_fill_manual(values =my_palette_large[1:10])+
  theme_classic() +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 14),#设定分页图的标题属性
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 10),
        legend.text = element_text(face = 'bold',color = 'black',size = 10),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1)) +
  labs(x='Position', y='Coverage') +
  ggtitle("Phi6_S")

graph2ppt(Phi6, file="c2",width = 9, aspectr=sqrt(2), append =TRUE)

############## FigS1 Qubit vs Nanodrop ################

concentration <- read.csv("concentration_qubit_nanodrop.csv",sep=',',header = T)
concentration1 <- concentration[,-4]
concentration1 <- melt(concentration1)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

concentration2 <- data_summary(concentration1,varname='value',groupnames = c("phage","variable"))

my_comparisons <- list(c("T4","T4"))

figs1 <- concentration2 %>% ggplot(aes(x=phage, y= value, fill = variable)) +  
  geom_bar(stat = "identity", position = "dodge", width = .5) +
  scale_colour_manual(values =my_fun_col_palettes[4:5])+
  scale_fill_manual(values =my_fun_col_palettes[4:5])+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.5)) +
  #facet_wrap(~Category, nrow =1, scales = "free") +
  theme_classic() +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 5),
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 5),
        legend.text = element_text(face = 'bold',color = 'black',size = 5),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1)) +
        #stat_compare_means(value ~ phage, method = "t.test") +
ggtitle("Original_phage")

graph2ppt(figs1, file="Fig4A",paper= "A4",width = 15, height = 5,vector.graphic = TRUE)



################################ Pearson correlation analysis between expected and observed (%)
#figS2A <- ggplot(TPM,aes(x=Ex,y=real)) +
#  geom_point(aes(color = Phage, shape = Phage)) +
#  geom_smooth(aes(color = Phage), method = "lm", se = FALSE, fullrange = TRUE) +
#  #facet_wrap(~M) +
#  scale_color_manual(values = my_fun_col_palettes) +
#  scale_fill_manual(values =  my_fun_col_palettes) +
#  ggpubr::stat_cor(aes(color = Phage), label.x = 7.5, method = "pearson", digits = 4) +
#  theme_classic() +
#  theme(strip.text = element_text(face = 'plain',color = 'black',size = 14),
#        axis.title.y = element_text(face = 'plain',color = 'black',size = 14),
#        axis.title.x = element_text(face = 'plain',color = 'black',size = 14),
#        axis.text.y = element_text(face = 'plain',color = 'black',size = 14),
#        axis.text.x = element_text(face = 'plain',color = 'black',size = 14),
#        legend.text = element_text(face = 'plain',color = 'black',size = 14),
#        legend.position="right",
#        panel.grid = element_blank(),
#        legend.background = element_blank(),
#        strip.background = element_blank(),
#        axis.text.x.bottom = element_text(angle = 0,hjust = 1, vjust = 1)) +
#  ggtitle("Mock_C with T4")
#graph2ppt(figS2A, file="./figS2_Mock_C with T4",paper= "A4",width = 5, height = 3,scaling = 100,vector.graphic = TRUE)

library("ggpubr")
#TPM <- read.csv("C:/Users/nwb188/OneDrive - University of Copenhagen/PhD/SSL/SSL_data/R/Fig2S_DNA_NXT077_pearsoncor-T4_new_dsDNA.csv",sep = ',', header = T)
my_cor <- c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#91D1C2FF", "#7E6148FF")

TPM2 <- TPM %>% filter(!Method %in% c("ExpectedA","ExpectedB","ExpectedC"))
TPM$Library.method = factor(TPM$Library.method, levels=c('SSLR', 'xGen', 'MDA_0.5h','MDA_1.5h', 'Nextera'))
Fig2C <- TPM %>% filter(T4=="Yes") %>%filter(!Method %in% c("SSLR100","ExpectedA","ExpectedB","ExpectedC")) %>%
         filter(!Phage %in% c("M13mp18","PhiX174")) %>%
          ggscatter(x="Ex",y="real",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Expected (%)", ylab = "Observed (%)") +
  geom_point(aes(color = Phage),size = 3) +
  #geom_smooth(aes(color = Phage), method = "lm", se = FALSE, fullrange = TRUE) +
  facet_wrap(~Method,scales = "free_x", ncol = 6) +
  scale_color_manual(values = my_palette_large[1:10]) +
  scale_fill_manual(values =  my_palette_large[1:10]) +
  #ggpubr::stat_cor(aes(color = Phage), label.x = 7.5, method = "pearson", digits = 4) +
  theme_classic() +
  theme(strip.text = element_text(face = 'plain',color = 'black',size = 14),
        axis.title.y = element_text(face = 'plain',color = 'black',size = 14),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 14),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 14),
        axis.text.x = element_text(face = 'plain',color = 'black',size = 14),
        legend.text = element_text(face = 'plain',color = 'black',size = 14),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 0,hjust = 1, vjust = 1)) +
  ggtitle("dsDNA Mock_pearson with T4")

graph2ppt(Fig2C, file="./FV_result/dsDNA Mock_pearson with T4",paper= "A4",width = 15, height = 5,scaling = 100)

################################ Pearson correlation analysis between expected and observed (%)
TPM3 <- TPM %>% filter(!Method %in% c("ExpectedA","ExpectedB","ExpectedC"))
TPM$Library.method = factor(TPM$Library.method, levels=c('SSLR', 'xGen', 'MDA_0.5h','MDA_1.5h', 'Nextera'))
Fig2C <- TPM %>% filter(T4=="Yes") %>%filter(!Method %in% c("SSLR100","ExpectedA","ExpectedB","ExpectedC")) %>%
  filter(!Phage %in% c("M13mp18","PhiX174")) %>%
  ggscatter(x="Ex",y="real",
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Expected (%)", ylab = "Observed (%)") +
  geom_point(aes(color = Phage),size = 3) +
  #geom_smooth(aes(color = Phage), method = "lm", se = FALSE, fullrange = TRUE) +
  facet_wrap(~Method,scales = "free_x", ncol = 6) +
  scale_color_manual(values = my_palette_large[1:10]) +
  scale_fill_manual(values =  my_palette_large[1:10]) +
  #ggpubr::stat_cor(aes(color = Phage), label.x = 7.5, method = "pearson", digits = 4) +
  theme_classic() +
  theme(strip.text = element_text(face = 'plain',color = 'black',size = 14),
        axis.title.y = element_text(face = 'plain',color = 'black',size = 14),
        axis.title.x = element_text(face = 'plain',color = 'black',size = 14),
        axis.text.y = element_text(face = 'plain',color = 'black',size = 14),
        axis.text.x = element_text(face = 'plain',color = 'black',size = 14),
        legend.text = element_text(face = 'plain',color = 'black',size = 14),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 0,hjust = 1, vjust = 1)) +
  ggtitle("dsDNA Mock_pearson wit T4")

graph2ppt(Fig2C, file="./R/DNA_Mock/FV_result/dsDNA Mock_pearson with T4",paper= "A4",width = 15, height = 5,scaling = 100)
TPM4 <- TPM[-c(8,17,26,35,44,53,62,71,80,89,98,107,116,125,134,143,152,161,170,179,188,197,206,215,224,233,242,251,260,269,278,287,296,305,314,323),]
TPM4 <- TPM[-c(8,17,26,35,44,53,62,71,80,89,98,107,116,125,134,143,152,161,170,179,188,197,206,215,224,233,242,251,260,269,278,287,296,305,314,323,332,341,350,359,368,377),]
############## FigS4 MDA Qubit vs Nanodrop ################

concentration <- read.csv("MDA_concentration_qubit_nanodrop.csv",sep=',',header = T)
concentration1 <- melt(concentration,id.vars = ("phage"), measure.vars=c("Nanodrop", "Qubit"), variable.name="Concentration")
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
library("plyr")
library("dplyr")
concentration2 <- data_summary(concentration1,varname="value",groupnames = c("phage","Concentration"))

my_comparisons <- list(c("T4","T4"))

figs2 <- concentration2 %>% ggplot(aes(x=phage, y= value, fill = Concentration)) +  
  geom_bar(stat = "identity", position = "dodge", width = .5) +
  scale_colour_manual(values =my_fun_col_palettes)+
  scale_fill_manual(values =my_fun_col_palettes)+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.5)) +
  #facet_wrap(~Category, nrow =1, scales = "free") +
  theme_classic() +
  geom_text(aes(label = value, y = value-2*sd), position = position_dodge( 0.5),vjust=0) +
  theme(strip.text = element_text(face = 'bold',color = 'black',size = 5),
        axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 5),
        legend.text = element_text(face = 'bold',color = 'black',size = 5),
        legend.position="right",
        panel.grid = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1, vjust = 1)) +
  #stat_compare_means(value ~ phage, method = "t.test") +
  ggtitle("Original_phage")

graph2ppt(figs2, file="figs2",paper= "A4",width = 15, height = 5,vector.graphic = TRUE)



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
################
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
################
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
##############
data <- read.table(header=TRUE, text='
 Subject RoundMono SquareMono RoundColor SquareColor
                   1        41         40         41          37
                   2        57         56         56          53
                   3        52         53         53          50
                   4        49         47         47          47
                   5        47         48         48          47
                   6        37         34         35          36
                   7        47         50         47          46
                   8        41         40         38          40
                   9        48         47         49          45
                   10        37         35         36          35
                   11        32         31         31          33
                   12        47         42         42          42
                   ')

# Convert it to long format
library(reshape2)
data_long <- melt(data=data, id.var="Subject",
                  measure.vars=c("RoundMono", "SquareMono", "RoundColor", "SquareColor"),
                  variable.name="Condition")
names(data_long)[names(data_long)=="value"] <- "Time"

# Split Condition column into Shape and ColorScheme
data_long$Shape <- NA
data_long$Shape[grepl("^Round",  data_long$Condition)] <- "Round"
data_long$Shape[grepl("^Square", data_long$Condition)] <- "Square"
data_long$Shape <- factor(data_long$Shape)

data_long$ColorScheme <- NA
data_long$ColorScheme[grepl("Mono$",  data_long$Condition)] <- "Monochromatic"
data_long$ColorScheme[grepl("Color$", data_long$Condition)] <- "Colored"
data_long$ColorScheme <- factor(data_long$ColorScheme, levels=c("Monochromatic","Colored"))

# Remove the Condition column now
data_long$Condition <- NULL


datac <- summarySEwithin(data_long, measurevar="Time", withinvars=c("Shape","ColorScheme"), idvar="Subject")
datac
#>    Shape   ColorScheme  N     Time Time_norm       sd        se        ci
#> 1  Round       Colored 12 43.58333  43.58333 1.212311 0.3499639 0.7702654
#> 2  Round Monochromatic 12 44.58333  44.58333 1.331438 0.3843531 0.8459554
#> 3 Square       Colored 12 42.58333  42.58333 1.461630 0.4219364 0.9286757
#> 4 Square Monochromatic 12 43.58333  43.58333 1.261312 0.3641095 0.8013997

library(ggplot2)
ggplot(datac, aes(x=Shape, y=Time, fill=ColorScheme)) +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=Time-ci, ymax=Time+ci)) +
  coord_cartesian(ylim=c(40,46)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_y_continuous(breaks=seq(1:100)) +
  theme_bw() +
  geom_hline(yintercept=38) 



#################gene comparison with different library preparations#################

install.packages("CMplot", lib = ("C:/Program Files/R/R-4.1.3/library"))

library("CMplot")

SNP <- read.delim("iSeq007_4_S4_L001.txt",sep="\t",header=T)



#data(pig60K)   #calculated p-values by MLM
#data(cattle50K)   #calculated SNP effects by rrblup
#head(pig60K)
SNP$MDA_1.5h_coverage <- SNP$MDA_1.5h_coverage + 1.1
SNP$MDA_1.5h_coverage <- log10(SNP$MDA_1.5h_coverage)/4
#SNP[is.na(SNP)] <- 0
SNP[is.na(SNP)] <- 0.00001
write.table(SNP,file = "SNP1.txt",sep="\t")

p5 <- CMplot(SNP,type="p",plot.type="c",r=0.4,cir.legend=TRUE,LOG10=FALSE,threshold=NULL,
         outward=FALSE,cir.legend.col="black",cir.chr.h=1.0,chr.den.col="black",
         memo="",dpi=300,file.output= FALSE,verbose=TRUE,width=10,height=10)
graph2tif(p5, file="./SSL/SSL data/R/test", width = 10, height = 10, scaling = 100)
graph2ppt(p5, file="./Fig3A_coverage")

CMplot(SNP,plot.type="q",col=matrix(c("#4DAF4A","dodgerblue4","deepskyblue","dodgerblue1", 
                                      "olivedrab3", "darkgoldenrod1")),threshold=1e-6,
       ylab.pos=2,signal.pch=c(19,6,4),signal.cex=1.2,signal.col="red",conf.int=TRUE,box=FALSE,multracks=
         TRUE,cex.axis=2,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,ylim=c(0,8),width=5,height=5)

LOG10=FALSE,threshold=NULL,
