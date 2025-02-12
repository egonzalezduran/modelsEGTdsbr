####################### R Code for statistical models (GLM) ##############################################
"Suppression of plastid-to-nucleus gene transfer by DNA double-strand break repair "                
Gonzalez-Duran, E and Bock, R (2025)                                                               
Version 06.02.25 by Enrique Gonzalez-Duran                                                         
R version 3.5.3 accessed through R studio  2023.12.0 Build 369                                      
Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                            
                                                                                                       
#####################################################################################################   

This code is set to work in R version 3.5.3.

List of necessary additional packages:

MASS        v. 7.3-51.1
dplyr       v. 0.8.0.1
ggplot2     v. 3.1.1
ciTools     v. 0.5.1
insight     v. 0.2.0
pscl        v. 1.5.2
MuMIN       v. 1.43.6
rstudioapi  v. 0.10

In order to run the code, just put all the files in a single folder of your choosing and run the Models for DNA repair effect over EGT.R using R studio.
The program should be able to recognize the path where the R script is, set it as the working directory, and use the datafiles from there. 
Output files will be generated in the working directory.

