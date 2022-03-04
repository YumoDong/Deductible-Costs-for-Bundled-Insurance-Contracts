# Deductible-Costs-for-Bundled-Insurance-Contracts

This directory provides data and codes to reproduce the results in the following paper:

Dong, Y., Frees, E., Huang, F.. (2021). Deductible Costs for Bundled Insurance Contracts, Working paper.

Software Dependency
===================
The plots/tables were produced under the following software
1. RStudio 1.4.1103
2. Operating system: Windows 10

[Important Notes]:

1. The deductible cost in this paper is calculated by the Monte Carlo simulation method. The time of the simulation process is affected by the number of simulation, so we provide the source codes of the Monte Carlo simulation process (codeforApp1.R, codeforApp2.R, codeforApp3.R, and codeforsimulationstudy.R) and also the simulated data sets (simulation1.RData, simulation2.RData, simulation3.RData, and simulationstudy.RData).
2. Deductible-Costs-YMD_blinded.rmd is the RMarkdown file that generates the main body of the paper including all the figures and tables. Thus, to reproduce any result shown in the paper, please refer to this file for the corresponding R code.
3. Please read the following Guidlines for Reproduction to see the details of reproducing the result in each section of the paper.

Guidlines for Reproduction
=========================
In Deductible-Costs-YMD_blinded.rmd, we use the simulated data sets directly to conduct analysis. If you are NOT interested in how to simulate data sets, then you only need to load the data sets and run the code in Deductible-Costs-YMD_blinded.rmd for reproduction. If you are interested in how to get those simulated data sets, please read the following code of simulation proccess.

#################################################################

Application 1: Wisconsin Property Fund (Commercial)

"codeforApp1.R" is the source code for "simulation1.RData". Note we directly use the results from Jed's paper and save it in "outsamplepolicyholders.RData". If you are interested in those results, click this link for the data, code, and the paper: https://sites.google.com/a/wisc.edu/jed-frees/. Then check the section called "Wisconsin Property Fund".

##################################################################

Application 2: Insurance Company Loss and Expense (Reinsurance)

"codeforApp2.R" is the source code for "simulation2.RData".

##################################################################

Application 3: Quotes of a Pet Insurance Contract (Personal)

"codeforApp3.R" is the source code for "simulation3.RData".

##################################################################

SIMULATION STUDY

"codeforsimulationstudy.R" is the source code for "simulationstudy.RData".
