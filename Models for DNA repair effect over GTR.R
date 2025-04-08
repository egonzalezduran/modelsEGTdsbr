####################### R Code for statistical models (GLM) ##############################################
## "Suppression of plastid-to-nucleus gene transfer by DNA double-strand break repair "
## Gonzalez-Duran, E., Kroop, X., Schadach, A. and Bock, R. (2025)
## Version 08.04.25 by Enrique Gonzalez-Duran
## R version 4.3.3 
## Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                           ###  
## PART 1                                                                                              ###   
###################################################################################################### ###   

##This code is set to work in R version 4.3.3.

##List of necessary additional packages:

# MASS        v. 7.3-60.0.1
# dplyr       v. 1.1.4
# ggplot2     v. 3.5.1
# ciTools     v. 0.6.1
# insight     v. 0.20.0
# rstudioapi  v. 0.16.0
# multcomp    v. 1.4-25

wants <- c("MASS","dplyr","ggplot2","insight","ciTools","rstudioapi","multcomp") ## searches for and installs the packages
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

library(MASS) 
library(dplyr) 
library(ggplot2)
library(ciTools)
library(insight)
library(rstudioapi)
library(multcomp)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets the location of the file as the working directory, works in R Studio
getwd()
################### Model 1 ########################

FGT1table <- read.table("SomaticFGT1.txt", header = TRUE, na.strings="z")
FGT1df <- as.data.frame(FGT1table)
FGT1df

M1 <-glm(FGT ~ 0 + Genotype + SELt:Genotype, weights= sqrt(TOTev), family=gaussian, data = FGT1df, contrasts= list(Genotype=contr.treatment(c("lig41p","polqpol1","pRB98"),base=3)))
confint(M1)
summary(M1)

FGT1wci<-add_ci(FGT1df,M1,alpha=0.05)
FGT1wci

coefsM1<- as.vector(coef(M1))

coefFGT1df<-data.frame(matrix(data=NA, nrow=6, ncol= 0),row.names = c("Genotypelig4-1","GenotypePolqpol-1","GenotypepRB98","Slope lig4-1p","Slope polqpol1","Slope pRB98"), stringsAsFactors = TRUE)
coefFGT1df$coefficientsFGT1<- as.vector(coef(M1)/5000000)
coefFGT1df

write.table(coefFGT1df, file= "M1 coefficients.txt", sep= "\t", dec = ".", row.names = TRUE, col.names = TRUE)
sink("SummaryM1.txt")
print(summary(M1))
sink()  # returns output to the console


scientific1 <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text=l)
}
plot.M1 <- ggplot(data= FGT1wci, aes(x = SELt, y = FGT, group = Genotype, colour=Genotype)) +
  scale_x_continuous(name="Time on selection (days)", limits=c(0, 140),breaks=c(0,20,40,60,80,100,120))+
  scale_y_continuous(name="FGT (events per cell)", limits=(c(-2.5,15)/5000000), breaks = (c(0 , 2.5,5,7.5,10,12.5,15)/5000000),labels=scientific1 ) +
  geom_line(aes(x= SELt, y=pred/5000000, group=Genotype, colour=Genotype),linetype = "longdash",alpha=0.8,size=1 )+
  geom_point(aes(y= FGT/5000000,fill=Genotype), alpha = 0.8, size = 2, shape=21) +
  geom_line(data=FGT1wci, aes(x= SELt, y=pred/5000000, group=Genotype, colour=Genotype),alpha=0.6,size=1)+
  geom_ribbon(data=FGT1wci,aes(ymax=UCB0.975/5000000, ymin= LCB0.025/5000000, group= Genotype, fill=Genotype),colour=NA, alpha=0.4) +
  geom_point(aes(y= 1/5000000,x=60,fill='Stegemann et al. 2003'), alpha = 1, size = 2, shape=21) +
  scale_fill_manual(values=c("#8400bd","#faaa44","#1F78B4","#AAFF00")) +
  scale_colour_manual(values=c("#8400bd","#faaa44","#1F78B4","#AAFF00")) +
  theme_bw()  
plot.M1

SM1<-summary(M1)
PM1 <- get_parameters(M1)
CIM1 <- confint.default(M1)
PM1[,3]<-SM1$coefficients[,2]
PM1 <- cbind(PM1,CIM1)
PM1[,6]<-SM1$coefficients[,3]
PM1[,7]<-SM1$coefficients[,4]

Master.table.PM.hlp<-as.data.frame(bind_rows(PM1))
colnames(Master.table.PM.hlp) <- c("Parameter","Estimate","Standard.Error","CI95.lower.bound","CI95.upper.bound", "Wald.z-statistic", "P-value")
Master.table.PM.hlp[,1] <- c("M1 Incercept lig4-1p","M1 Incercept polqpol1", "M1 Intercept Nt-RB98", "M1 Slope lig4-1p", "M1 Slope polqpol-1", "M1 Slope Nt-RB98")
Master.table.PM.hlp

Master.table.PM.cell <-Master.table.PM.hlp
Master.table.PM.cell[,2] <-Master.table.PM.hlp[,2]/5000000
Master.table.PM.cell[,3] <-Master.table.PM.hlp[,3]/5000000
Master.table.PM.cell[,4] <-Master.table.PM.hlp[,4]/5000000
Master.table.PM.cell[,5] <-Master.table.PM.hlp[,5]/5000000
colnames(Master.table.PM.cell) <- c("Parameter","Estimate","Standard.Error","CI95.lower.bound","CI95.upper.bound", "Wald.z-statistic", "P-value")
Master.table.PM.cell

sink(file= "M1 FGT per cell over time.txt")
Master.table.PM.cell
sink()
write.table(Master.table.PM.cell, file="M1 FGT per cell over time.csv", sep = "  ")

################### Model 2 ########################

FGT2table <- read.table("SomaticFGT2.txt", header = TRUE, na.strings="z")
FGT2df <- as.data.frame(FGT2table)
FGT2df

M2 <-glm(FGT ~ 0 + Genotype + SELt:Genotype , weights= sqrt(TOTev), data = FGT2df, family = gaussian, contrasts= list(Genotype=contr.treatment(c("lig42","polqpol2","pRB98"),base=3)))
confint(M2)
summary(M2)
FGT2wci<-add_ci(FGT2df,M2,alpha=0.05)
FGT2wci

coefsM2<- as.vector(coef(M2))

coefFGT2df<-data.frame(matrix(data=NA, nrow=6, ncol= 0),row.names = c("Genotypelig4-2","GenotypePolqpol-2","GenotypepRB98","Slope lig4-1","Slope polqpol1","Slope pRB98"), stringsAsFactors = TRUE)
coefFGT2df$coefficientsFGT2<- as.vector(coef(M2)/5000000)
coefFGT2df

write.table(coefFGT2df, file= "M2 coefficients.txt", sep= "\t", dec = ".", row.names = TRUE, col.names = TRUE)
sink("SummaryM2.txt")
print(summary(M2))
sink()  # returns output to the console

plot.M2 <- ggplot(data= FGT2wci, aes(x = SELt, y = FGT, group = Genotype, colour=Genotype)) +
  scale_x_continuous(name="Time on selection (days)", limits=c(0, 140),breaks=c(0,20,40,60,80,100,120))+
  scale_y_continuous(name="FGT (events per cell)", limits=(c(0,15)/5000000), breaks = (c(0 , 2.5,5,7.5,10,12.5,15)/5000000),labels=scientific1 ) +
  geom_line(aes(x= SELt, y=pred/5000000, group=Genotype, colour=Genotype),linetype = "longdash",alpha=0.8,size=1 )+
  geom_point(aes(y= FGT/5000000,fill=Genotype), alpha = 0.8, size = 2, shape=21) +
  geom_line(data=FGT2wci, aes(x= SELt, y=pred/5000000, group=Genotype, colour=Genotype),alpha=0.6,size=1)+
  geom_ribbon(data=FGT2wci,aes(ymax=UCB0.975/5000000, ymin= LCB0.025/5000000, group= Genotype, fill=Genotype),colour=NA, alpha=0.4) +
  scale_fill_manual(values=c("#8400bd","#faaa44","#1F78B4")) +
  scale_colour_manual(values=c("#8400bd","#faaa44","#1F78B4")) +
  theme_bw()  
plot.M2


SM2<-summary(M2)
PM2 <- get_parameters(M2)
CIM2 <- confint.default(M2)
PM2[,3]<-SM2$coefficients[,2]
PM2 <- cbind(PM2,CIM2)
PM2[,6]<-SM2$coefficients[,3]
PM2[,7]<-SM2$coefficients[,4]

Master.table.PM2.hlp<-as.data.frame(bind_rows(PM2))
colnames(Master.table.PM2.hlp) <- c("Parameter","Estimate","Standard.Error","CI95.lower.bound","CI95.upper.bound", "Wald.z-statistic", "P-value")
Master.table.PM2.hlp[,1] <- c("M2 Incercept lig4-2","M2 Incercept polqpol2", "M2 Intercept Nt-RB98", "M2 Slope lig4-2", "Slope polqpol-2", "Slope Nt-RB98")
Master.table.PM2.hlp

Master.table.PM2.cell <-Master.table.PM2.hlp
Master.table.PM2.cell[,2] <-Master.table.PM2.hlp[,2]/5000000
Master.table.PM2.cell[,3] <-Master.table.PM2.hlp[,3]/5000000
Master.table.PM2.cell[,4] <-Master.table.PM2.hlp[,4]/5000000
Master.table.PM2.cell[,5] <-Master.table.PM2.hlp[,5]/5000000
colnames(Master.table.PM2.cell) <- c("Parameter","Estimate","Standard.Error","CI95.lower.bound","CI95.upper.bound", "Wald.z-statistic", "P-value")
Master.table.PM2.cell

sink(file= "M2 FGT per cell over time.txt")
Master.table.PM2.cell
sink()
write.table(Master.table.PM2.cell, file="M2 FGT per cell over time.csv", sep = "  ")

######### Model 3##############

GTpollen.df <-  as.data.frame(matrix(data=NA, nrow=7, ncol=5))
GTpollen.df[1,] <- c("RB98","none",3,154400,NA)
GTpollen.df[2,] <- c("lig4-1","lig4",12,89400,NA)
GTpollen.df[3,] <- c("lig4-2","lig4",12,72400,NA)
GTpollen.df[4,] <- c("polqpol-1","polq",9,94600,NA)
GTpollen.df[5,] <- c("polqpol-2","polq",11,108000,NA)
GTpollen.df[6,] <- c("polqcds-3","polq",25,65600,NA)
GTpollen.df[7,] <- c("polqhel-4","polq",22,62800,NA)
GTpollen.df[,3]<- as.numeric(GTpollen.df[,3])
GTpollen.df[,4]<- as.numeric(GTpollen.df[,4])
colnames(GTpollen.df) <- c("Genotype","Affected.gene.DSBR","GT.events","Total.seeds","FGT")
GTpollen.df$Genotype = factor(GTpollen.df$Genotype,
                                        levels=unique(GTpollen.df$Genotype))
GTpollen.df$Affected.gene.DSBR = factor(GTpollen.df$Affected.gene.DSBR,
                        levels=unique(GTpollen.df$Affected.gene.DSBR))


GTpollen.df[,5] <- GTpollen.df$GT.events / GTpollen.df$Total.seeds
GTpollen.df

comparison.vector1<-as.vector(as.factor(unique(GTpollen.df[,1])))
comparison.vector1


M3 <- glm(FGT ~ Genotype, weights = Total.seeds, family= binomial(link= "log"),contrast=list(Genotype=contr.treatment(comparison.vector1,base=1)),data=GTpollen.df)
summary(M3)

hpsetM3<-cftest(glht(M3), c("Genotypelig4-1","Genotypelig4-2","Genotypepolqpol-1","Genotypepolqpol-2","Genotypepolqcds-3","Genotypepolqhel-4"), test= Chisqtest())
hptestedcM3<-summary(hpsetM3, test=adjusted(type="holm"))
confint(glht(M3),c("Genotypelig4-1","Genotypelig4-2","Genotypepolqpol-1","Genotypepolqpol-2","Genotypepolqcds-3","Genotypepolqhel-4"),level=0.95,calpha = adjusted_calpha())
hptestedcM3

sink(file= "Hypothesis testing M3.txt")
hptestedcM3
sink()

#"holm" is the Holm-Bonferroni method for Family Wise Error rate (FWER) control. It is more powerful than Bonferroni, but less than FDR. 

#### Confidence Intervals
GTpollenCI95M3 <- as.data.frame(add_ci(GTpollen.df,M3, alpha = 0.05))
GTpollenCI95M3
GTpollenCI95M3<- transform(GTpollenCI95M3, FGT = as.numeric(FGT), 
                         GT.events  = as.numeric(GT.events),
                         Total.seeds  = as.numeric(Total.seeds),
                         pred  = as.numeric(pred),
                         UCB0.975  = as.numeric(UCB0.975),
                         LCB0.025  = as.numeric(LCB0.025))


mycolorsp <- c("#417dd1","#8400bd","#faaa44")

brks = seq(0, 8, 0.5)
lbs <- ifelse(brks%%2 == 0,format(brks/10000, digits = 1, nsmall = 1),"")

plot.M3 <- ggplot(GTpollenCI95M3, aes(x=Genotype, y= FGT,color= Affected.gene.DSBR)) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_jitter(size=3, position = position_dodge2(width= 0.2) , alpha=0.5, show.legend= TRUE) + 
  geom_errorbar(aes(ymin=LCB0.025, ymax=UCB0.975),width=0,linetype=1, show.legend= FALSE, size=1) +  #plots the CI
  scale_color_manual(values= mycolorsp) +
  scale_fill_manual(values= mycolorsp) +
  theme(axis.line = element_line(colour = "black")) +
  xlab("Genotype") +
  ylab("FGT (EGT events/seedling)") +
  scale_y_continuous(limits = c(0, 0.0006), expand= c(0,0),breaks = brks/10000, labels = lbs)
plot.M3

##########Plots###########


pdf(file="M1plot.pdf", height = 5, width =  8)
plot.M1 #here goes the plot code###
dev.off()

pdf(file="M2plot.pdf", height = 5, width =  8)
plot.M2 #here goes the plot code###
dev.off()

pdf(file="M3plot.pdf", height = 5, width =  8)
plot.M3 #here goes the plot code###
dev.off()


