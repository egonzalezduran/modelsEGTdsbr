# modelsEGTdsbr
####################### R Code for statistical models (GLM) ############################################## 

"Suppression of plastid-to-nucleus gene transfer by DNA double-strand break repair "
Gonzalez-Duran, E., Kroop, X., Schadach, A. and Bock, R. (2025)
Version 08.04.25 by Enrique Gonzalez-Duran
R versions 4.3.3 and 3.6.3 
Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany

##########################################################################################################

This code is split in two parts:

Models for DNA repair effect over GTR.R (part 1)
Slope comparisons.R (part 2)

Both scripts need to reed data from two data files:
SomaticFGT1.txt 
SomaticFGT2.txt

The part 1 contains most of the code and statistical tests. It works in R v.4.3.3 using the following packages:

# MASS        v. 7.3-60.0.1
# dplyr       v. 1.1.4
# ggplot2     v. 3.5.1
# ciTools     v. 0.6.1
# insight     v. 0.20.0
# rstudioapi  v. 0.16.0
# multcomp    v. 1.4-25

In order to run the code, just put all the files in a single folder of your choosing and run the Models for DNA repair effect over EGT.R using R studio. The program should be able to recognize the path where the R script is, set it as the working directory, and use the datafiles from there. Necessary packages are installed automatically. Output files will be generated in the working directory.

Part 2 needs to be run in R v.3.6.3. Part 2 requires R the phia package v.0.2-1, which does NOT work in R v.4.3.3 (The Testinteraction() function doesn't work as intended). 
To run the code, place the slope comparison.R file together with the SomaticFGT1.txt and SomaticFGT2.txt files in the working directory. Run the script, which installs the package in its correct version. Recover the output files from the working directory.
