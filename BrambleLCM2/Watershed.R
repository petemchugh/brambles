####################################################
# ISEMP Watershed Production Model
# R-Code Written by Matt Nahorniak, South-Fork Research
# matt@southforkresearch.org
# (541) 740-5487
# Pete McHugh's John Day 2015-2016 tweaks are included
# peter.a.mchugh@gmail.com
#####################################################
rm(list = ls(all=TRUE)) #Clean up shop before beginning
#setwd("C:/Users/etal/Desktop/Pete/S1 (msS3a)")
#install.packages("MCMCpack")
#install.packages("VGAM")
#install.packages("rprojroot")
library(MCMCpack, quietly=T) # for Dirichlet (multivariable beta) distributions
library(VGAM, quietly=T) # Library with posnorm functions (positiive only normal)


# Load R-Scripts in Different Files
#source("Watershed_ReadData.r")
source('~/BrambleCloudLCM2/Watershed_ReadData.R', encoding = 'UTF-8')
#source("Watershed_MonteCarlo.r")
source('~/BrambleCloudLCM2/Watershed_MonteCarlo.R', encoding = 'UTF-8')
#source("Watershed_BevHolt.r") # Version for other applications
#source("Watershed_BevHolt_NFLewisMAupd.r") # Version with supplementation processing activated for NF Lewis
source('~/BrambleCloudLCM2/Watershed_BevHolt_NFLewisMAupd.R', encoding = 'UTF-8')
#source("Watershed_Post_Processing.r")
source('~/BrambleCloudLCM2/Watershed_Post_Processing.R', encoding = 'UTF-8')
#source("Watershed_JohnDayOutputs.r") #Only use this if you're interested; it's for 1 pop and hasn't been vetted for Chinook
source('~/BrambleCloudLCM2/Watershed_JohnDayOutputs.R', encoding = 'UTF-8')

# Read Header and Input Files
print("Reading Header Files")
header<- Read.Header("Watershed_Header_File.csv") 
Inputs<-Read.Input.File(header)

# Initialze some Global Variables
# NR is global variable containing all counts N in run R
#NR = Abundance for site k, stage i, time t, life stage g, run R
#NR_F = fraction of NR that are females
#N5R = NR except just for pre-smolts
#N5R_F = fraction of N5R that are females (a fraction)
#SmoltsR = pre-smolts attempting to smolt (# pre Bev-Holt survival)
NR=with(header, {array(rep(0,K*I*Tr*G*R),c(K,I,Tr,G,R))})
NR_F=with(header, {array(rep(0,K*I*Tr*G*R),c(K,I,Tr,G,R))}) 
N5R=with(header, {array(rep(0,K*I5*Tr*G*R), c(K,I5,Tr,G,R))})
N5R_F=with(header, {array(rep(0,K*I5*Tr*G*R), c(K,I5,Tr,G,R))})
Candidate_Smolt_NR = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})      
Female_Spawners_NR = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
###***Pete May 2015 Array Additions***
Male_Spawners_NR = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
AdultReturnsSupp = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
Bonneville_NT = with(header, {array(rep(0, K*10*Tr*G*R), c(K,10,Tr,G,R))})
Candidate_Smolt_ByAge = with(header, {array(rep(0,K*I5*Tr*G*R), c(K,I5,Tr,G,R))})
Resident_Spawners_NR=with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
ResidentN=with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
ResidentSexRatio=with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
Prod=with(header, {array(rep(0,K*I*Tr*R),c(K,I,Tr,R))})
Cap=with(header, {array(rep(0,K*I*Tr*R),c(K,I,Tr,R))})
outstuff = rep(NA,79) #goes with CreateJohnDayOut()
###***Pete May 2015 Array Additions***
Spawners_NR = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
Escapement_NR = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})
Female_Escapement_NR = with(header, {array(rep(0, K*Tr*G*R), c(K,Tr,G,R))})


####################################################################
# Here we loop through the code one for each monte-carlo simulation
# First we set the variables, then we add stochasticity via a monte-
# carlo code, then run the Beverton-Holt algorithm
r=1
for (r in 1:header$R) {
  cat(paste('Iteration #',r), "\n")
  variables=Initialize.Variables(header)
  parameters=MonteCarlo(Inputs, variables, header)
  results= BevHolt(parameters, variables, header) 
  
  # Store results for Each Iteration
  NR[,,,,r]=results$N
  N5R[,,,,r]=results$N5
  Female_Spawners_NR[,,,r] = results$Female_Spawners
  Spawners_NR[,,,r] = results$Spawners
  Escapement_NR[,,,r] = results$Escapement
  Female_Escapement_NR[,,,r] = results$Female_Escapement
  Candidate_Smolt_NR[,,,r] = results$Candidate_Smolt
  AdultReturnsSupp[,,,r]=results$AdultReturns
  #Pass a bunch of stuff to get summarized by CreateJohnDayOut -- Pete July 2015 Addition
  #Comment the next two lines out if you don't want this...
  #WARNING--this piece was developed for a single site, natural-origin fish
  #   only situation; thus it shouldn't be used for other systems
  outstuff<-CreateJohnDayOut(r,results,header,parameters,outstuff)
  
} # End of Loop through MC Iteration
################################################################

#WriteJohnDayOut(outstuff) #Pete July 2015 Addition
#gsub<-'all' #or anything you'd like, e.g., gsub<-c("Natural","H1")
gsub<-c('all')
gtypeIn<-c("Natural","H1")
WriteJohnDayOut(outstuff,'sub',gtypeIn,FALSE,21,60,'Base')#Pete July 2016 version
#WriteJohnDayOut <- function(outstuff,which,nat,trim,first,last,S) 
#which = 'all'/'sub'; this isolates it to just 'johnday out' vs all the N5, N arrays
#nat = TRUE/FALSE; this limits to just Gtype = Natural
#trim = TRUE/FALSE;  this allows burnin or any other slop to be trimmed here
#first & last = first/last years if trim = TRUE; where to trim if trim = TRUE 
#S is character string name for the scenario, e.g., S1, S2, S3, whatever you like...
# WARNING--this piece was developed for a single site, natural-origin fish
#   only situation; thus it shouldn't be used for other systems
# filtdat<-read.table(file.choose(),header==TRUE,sep=",")
# ddd<-subset(filtdat,Gtype=="Natural")
# ddd<-subset(filtdat,yr>20)
# ddd<-subset(filtdat,yr<=60)
# write.csv(ddd,"NatRerunS4_tryagain.csv")

#####################################################
# !!!! Code below was commented out to faciliate bramblecloud runs...
# PostProcessing(header, NR, N5R, Female_Spawners_NR) 
# 
# ####################################################################
# # Generate plot for coefficient of variation vs number of MC iterations
# 
# if (header$R > 1) {
#   if (header$K==1) {mar = 2} else {mar = 3}
#   r.range=2:200
#   escapement_Mean = mean(apply(Escapement_NR[,header$Tr-1,,],mar,sum)) #to correct error, Pete changed apply argument 2 to 2 (from 3)
#   escapement_SD = sd(apply(Escapement_NR[,header$Tr-1,,],mar,sum))
#   escapement_se_by_r = escapement_SD/sqrt(r.range)
#   escapement_COV_by_r = escapement_se_by_r/ escapement_Mean
#   presmolts_Cand_mean = mean(apply(Candidate_Smolt_NR[,header$Tr-1,,],mar,sum))
#   presmolts_Cand_SD = sd(apply(Candidate_Smolt_NR[,header$Tr-1,,],mar,sum))
#   presmolts_Cand_se_by_r= presmolts_Cand_SD/sqrt(r.range)
#   presmolts_Cand_COV_by_r= presmolts_Cand_se_by_r /  presmolts_Cand_mean
#   
#   plot(r.range, escapement_COV_by_r, type="l", 
#   main="Coefficient of Variation for Total Escapement and
#   Total Pre-Smolts, by Number of MC Iterations" ,
#   xlab= "Number of MC Iterations",
#   ylab="Coefficient of Variation",
#   ylim=c(0, 1))
#   lines(r.range, presmolts_Cand_COV_by_r, lt=2)
#   legend("topright", c("Escapement", "Smolts (Pre-Dam Passage)"), lt=c(1,2))
# }

# results$Rainbow_N[1,,1:25,1]
# results$Rainbow_N_F[1,,1:25,1]








