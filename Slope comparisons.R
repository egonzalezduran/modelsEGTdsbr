####################### R Code for statistical models (GLM) ##############################################
## "Suppression of plastid-to-nucleus gene transfer by DNA double-strand break repair "
## Gonzalez-Duran, E., Kroop, X., Schadach, A. and Bock, R. (2025)
## Version 08.04.25 by Enrique Gonzalez-Duran
## R version 3.6.3 
## Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                           ###  
## PART 2                                                                                              ###   
###################################################################################################### ###   


### Part 2
### This piece of code concerns only the comparison of the GTR (slopes of genotypes over time of selection) which are featured in Supplementary Table 1
### This comparison depends on the testInteractions function of the Phia package, which does NOT work in R.4.3.3. 
### To obtain the pvalues of Supplementary Table 1, evaluate the Chisq statistic in another program (e.g. using the =CHISQ.DIST.RT function in Microsoft Excel) 

install.packages("https://cran.r-project.org/src/contrib/Archive/phia/phia_0.2-1.tar.gz", repos = NULL, type = "source")
library(phia)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets the location of the file as the working directory, works in R Studio
getwd()

FGT1table <- read.table("SomaticFGT1.txt", header = TRUE, na.strings="z")
FGT1df <- as.data.frame(FGT1table)
FGT1df

FGT2table <- read.table("SomaticFGT2.txt", header = TRUE, na.strings="z")
FGT2df <- as.data.frame(FGT2table)
FGT2df

M1 <-glm(FGT ~ 0 + Genotype + SELt:Genotype, weights= sqrt(TOTev), family=gaussian, data = FGT1df, contrasts= list(Genotype=contr.treatment(c("lig41p","polqpol1","pRB98"),base=3)))
M2 <-glm(FGT ~ 0 + Genotype + SELt:Genotype , weights= sqrt(TOTev), data = FGT2df, family = gaussian, contrasts= list(Genotype=contr.treatment(c("lig42","polqpol2","pRB98"),base=3)))

sink(file= "GTR comparisons in M1 and M2.txt")
testInteractions(M1, slope="SELt")
testInteractions(M2, slope="SELt")
sink()

###For exact p-values, evaluate the Chisq value in Excel using =CHIDIST with df=1### 