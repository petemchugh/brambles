###################################################
# This File contains two functions:
#
# CreateJohnDayOut(X1,X2,X3,X4...) which computes summary statistics for each iteration
# WriteJohnDayOut(...) which writes the CSV file containing everything...
# 
# This function is very straightforward, it just does a bunch of post-processing
# It has been moved here for tidiness reasons
# NOTE: at this time this has been designed for a single population and is steelhead centric

# Updated Sept 19 2016 to produce spawner to spawner summaries
# and with some slight re-ordering and re-naming to make it 
# more user-friendly to others (for Chris J., Rich Z, etc.)

###################################################

CreateJohnDayOut<-function(r,results,header,parameters,outstuff)
  {
      #This is a bunch of output bundling for making sure the model properly reflects the input data
      #It's tailored to the John Day O mykiss case, i.e., 1-3 salt returns and smolt ages 1-3
      #Also, it is currently structured for testing purposes only, i.e., single rep runs (testing/validation)
      #First is an adult check, including returns to bonneville, age comp, sex ratio on spawning grounds, etc.
      #Note this only covers a life history with up to 3 ocean ages, plus repeat spawners
      
  #Pete Feb2016 multipop, anywhere kk appears
  #gg<-3 kk<-1
  Gname<-c("Natural","H1","N.H2","H2","N-N.H3","N-H3","N.H3-H3","N.H3-H2",
                 "H3-H2", "H3","N.H3-N.H3")
  Male_Spawners_NR[,,,r] = results$Male_Spawners #for computing sex ratio of spawners
  Female_Spawners_NR[,,,r] = results$Female_Spawners
  Bonneville_NT[,,,,r] = results$NT #for computing age comp of adult spawners
  Candidate_Smolt_ByAge[,,,,r] = results$CandSmoltsByAge
  Resident_Spawners_NR[,,,r] =results$RainbowSpawners
  ResidentSexRatio[,,,r] =results$RainbowFemSpawners/results$RainbowSpawners
  NR[,,,,r]=results$N
  N5R[,,,,r]=results$N5
  Prod[,,,r]=results$p
  Cap[,,,r]=results$c
  
  for(gg in 1:header$G)
  {
      for(kk in 1:header$K)
      {
        Gtype<-rep(Gname[gg],header$Tr)
        SimYr<-seq(1:header$Tr)
        REP<-rep(r,length(SimYr))
        SubPop<-rep(kk,length(SimYr))
        SubPopName<-rep(header$Site.Names[kk],length(SimYr))
        ###***End-Pete May 2015 Addition***
        adults<-t(Bonneville_NT[kk,1:4,,gg,r]) #Pete Feb2016 Lewis function
        tot<-adults[,1]+adults[,2]+adults[,3]+adults[,4]
        pct1<-adults[,1]/tot
        pct2<-adults[,2]/tot
        pct3<-(adults[,3]+adults[,4])/tot
        adults<-cbind(AdultReturnsSupp[kk,,gg,r],adults,pct1,pct2,pct3)
        sexrat<-cbind(Female_Spawners_NR[kk,,gg,r],Male_Spawners_NR[kk,,gg,r],Male_Spawners_NR[kk,,gg,r]+Female_Spawners_NR[kk,,gg,r],
                      Female_Spawners_NR[kk,,gg,r]/(Male_Spawners_NR[kk,,gg,r]+Female_Spawners_NR[kk,,gg,r]))
        adults<-cbind(REP,SimYr,adults,sexrat,(Male_Spawners_NR[kk,,gg,r]+Female_Spawners_NR[kk,,gg,r])/(tot))
        #Second, some smolt checks/calcs.
        smolty<-t(Candidate_Smolt_ByAge[kk,,,gg,r])[,1:3] #Pete Feb2016 Lewis function
        totsmSAR<-smolty[,1]+smolty[,2]+smolty[,3]
        smolty<-smolty/t(Prod[kk,,,r])[,c(5)] # Must divide by survival from trap to JDA to get it in Sm/Sp units at the trap...
        totsm<-smolty[,1]+smolty[,2]+smolty[,3]
        pct1sm<-smolty[,1]/totsm
        pct2sm<-smolty[,2]/totsm
        pct3sm<-smolty[,3]/totsm
        smolty<-cbind(smolty,totsm,pct1sm,pct2sm,pct3sm)
        #Third, some combined SAR and sm/sp calcs
        SAR<-rep(NA,header$Tr)
        sm.per.sp<-SAR
        for(ot in 2:(header$Tr-3))
        {
          SAR[ot]<-(adults[ot+1,3]+adults[ot+2,4]+adults[ot+3,5])/totsmSAR[ot]  # SAR is JDA to BON
          sm.per.sp[ot]<-sum(smolty[ot+1,1],smolty[ot+2,2],smolty[ot+3,3])/(sexrat[ot,1]+sexrat[ot,2]) 
        }
        adults<-cbind(Gtype,SubPop,SubPopName,adults)
        stages<-t(NR[kk,,,gg,r])
        stages<-stages[,1:11]
        stagesN5<-t(N5R[kk,,,gg,r])
        stagesN5<-stagesN5[,1:3]
        rainbowdat<-cbind(Resident_Spawners_NR[kk,,gg,r],ResidentSexRatio[kk,,gg,r],
                          t(results$Rainbow_N[kk,2:3,,gg]))
        
        AdRiverSurv<-as.matrix(1-parameters$harvest.wild[kk,])
        AdOcnSurv<-t(parameters$Sr[kk,,])[,8:11]
        SpostSp<-t(parameters$Post_Spawn_Survival_Anadromous_M[kk,,,gg])[,1:3]
        Surv<-cbind(t(Prod[kk,,,r])[,c(2:3)],t(parameters$SR5[kk,,])[,1:3],t(Prod[kk,,,r])[,c(5)],AdOcnSurv,AdRiverSurv,SpostSp)
        Capac<-cbind(t(Cap[kk,,,r])[,1:5],t(parameters$C_ocean[kk,,])[,1])
        SmoltProbF<-t(parameters$N5.Psmolt_F[kk,,])[,1:3]
        SmoltProbM<-t(parameters$N5.Psmolt_M[kk,,])[,1:3]
        Mat<-t(parameters$Mat8Plus_F[kk,,])[,1:4]
        Fec<-t(parameters$Female_Fecundity[kk,,,gg])[,1:4]
        DemParm<-cbind(Surv,SmoltProbF,SmoltProbM,Mat,Fec,Capac)
        
        if(r==1&&kk==1&&gg==1)
        {outstuff<-cbind(adults,smolty,SAR,sm.per.sp,stages,stagesN5,rainbowdat,DemParm)} else
        {outstuff<-rbind(outstuff,cbind(adults,smolty,SAR,sm.per.sp,stages,stagesN5,rainbowdat,DemParm))}
        
      }
      
  }
      return(outstuff)
  }

#This has been moved to this program for tidiness reasons.
WriteJohnDayOut <- function(outstuff,which,gtypeIn,trim,first,last,S) 
  #which = 'all'/'sub'; gtypeIn = string, 'all' = everything, or c("Natural","H1", ... [see below] );
  #options for gtypeIn are gtypeIn<-c("Natural","H1","N.H2","H2","N-N.H3","N-H3","N.H3-H3","N.H3-H2","H3-H2", "H3","N.H3-N.H3")
  #trim = TRUE/FALSE; first & last = first/last years if trim = TRUE 
  
  ########################################################################
  # Part 1 -- the main 'JohnDayOutputs' file...
  {
      namescol<-c("Gtype","SubPopN","Name","rep","yr","AdultRet","OA1pass","OA2pass","OA3pass","OA4pass","OA1.pct",
                  "OA2.pct","OA3plus.pct","FemEsc","MaleEsc","TotEsc","pct.FemEsc","Convert")
      namescolsmolt<-c("1smolt","2smolt","3smolt","totsmolt","1sm.pct","2sm.pct","3sm.pct")
      namestages<-c("PassTot","Negg","Nfry","Nparr","Npsm","Nsm","Nad0","Nad1","Nad2","Nad3","Nad4","Cand1","Cand2","Cand3")
      nameDem<-c("Segg","Sfry","Sps1","Sps2","Sps3","Ssm","Sad1o","Sad2o","Sad3o","Sad4o","SadR","Spost1","Spost2","Spost3")
      nameDem<-c(nameDem,"Psm1f","Psm2f","Psm3f","Psm1m","Psm2m","Psm3m","Pad1","Pad2","Pad3","Pad4")
      nameDem<-c(nameDem,"Fec1","Fec2","Fec3","Fec4","Cegg","Cfy","CParr","Cpresm","Csmolt","Cocean")
      namesnames<-c(namescol,namescolsmolt,"SAR","sm.per.sp",namestages,"ResSp","pct.ResFem","ResEgg","ResFry",nameDem)
      Fname<-paste("JohnDayOutputs_",S,"_",Sys.time(),".csv")
      Fname<-gsub("-","",Fname)
      Fname<-gsub(" ","_",Fname)
      Fname<-gsub(":","",Fname)
      #length(namesnames)
	outstuff1<-as.data.frame(outstuff,stringsAsFactors=FALSE)
	sadd<-rep(S,dim(outstuff1)[1])
	outstuff1<-cbind(sadd,outstuff1)
	names(outstuff1)<-c("Scenario",namesnames)
	#str(outstuff1); head(outstuff1)
	#trim<-TRUE; nat<-TRUE; first <- 5; last <- 10; S <- "S1"
      if(length(gtypeIn)!=1|gtypeIn[1]!= 'all') {
          outstuff1 <- subset(outstuff1,Gtype %in% gtypeIn)
      } #Gtype <- c(rep("Natural",2),rep("H2",2)) gtypeIn = c("Natural","H1") G<-Gtype subset(G,Gtype %in% gtypeIn)
	#str(outstuff1)
      if (trim == TRUE) {
	      outstuff1$yr<-as.numeric(outstuff1$yr)
	      keep<-(outstuff1$yr>=first & outstuff1$yr<=last)		
        outstuff2 <- as.data.frame(outstuff1[keep,]) #str(outstuff2)
      } else
      {
        outstuff2<-outstuff1
      }
  
      #re-order to make this a more sensible layout, a penalty for poor planning on my part :(
  
      outstuff3<-outstuff2[,c(1:7,29,8:28,30:80)] #df2[,c(1,3,2,4)]
	    namesnames2<-c(namesnames[c(1:6)],namesnames[c(28)],namesnames[c(7:27)],namesnames[c(29:79)])
      #names(outstuff) <- namesnames
      write.table(outstuff3,col.names=c("Scenario",namesnames2),
                  row.names=FALSE,sep = ",",file=Fname)
	

  
  
  ##############################################################################################
	##############################################################################################
	##############################################################################################
  #Additional code, added 20 Sept 2016, to estimate spawner/spawner ratios
  #comment out and use on an ad hoc basis if desired (just run on output from first part)
#   dat1<-read.table(file.choose(),header=TRUE,sep=",",stringsAsFactors=FALSE)
#   attach(dat1)
# 	#detach(dat1)
# 	outdim<-dim(dat1)[1]
#   minyr<-min(dat1$yr)
# 	maxyr<-max(dat1$yr)
#   kmax<-max(dat1$SubPopN)
#   gmon<-unique(dat1$Gtype)
#   scen<-unique(dat1$Scenario)
#   glen<-length(gmon)
#   maxage<-7 # this is total age, FW age + OA...you have to change this for Chinook & coho
# #new names
#   OA1sp<-TotEsc*OA1.pct #rep(NA,outdim)
# 	OA2sp<-TotEsc*OA2.pct #rep(NA,outdim)
# 	OA3sp<-TotEsc*OA3plus.pct #rep(NA,outdim)
# # #old names	
# #   OA1sp<-TotEsc*X1S.pct #rep(NA,outdim)1S.pct  2S.pct	3S+.pct
# # 	OA2sp<-TotEsc*X2S.pct #rep(NA,outdim)
# # 	OA3sp<-TotEsc*X3S..pct #rep(NA,outdim)
# 	BYsm1.N<-rep(NA,outdim)
# 	BYsm2.N<-rep(NA,outdim)
# 	BYsm3.N<-rep(NA,outdim)
# 	BYsm1.pct<-rep(NA,outdim)
# 	BYsm2.pct<-rep(NA,outdim)
# 	BYsm3.pct<-rep(NA,outdim)
# 	BYsm1.N<-rep(NA,outdim)
# 	BYsm2.N<-rep(NA,outdim)
# 	BYsm3.N<-rep(NA,outdim)
# 	BYsm1Ad1.N<-rep(NA,outdim)
# 	BYsm1Ad2.N<-rep(NA,outdim)
# 	BYsm1Ad3.N<-rep(NA,outdim)
# 	BYsm2Ad1.N<-rep(NA,outdim)
# 	BYsm2Ad2.N<-rep(NA,outdim)
# 	BYsm2Ad3.N<-rep(NA,outdim)
# 	BYsm3Ad1.N<-rep(NA,outdim)
# 	BYsm3Ad2.N<-rep(NA,outdim)
# 	BYsm3Ad3.N<-rep(NA,outdim)
# 	RecTA3<-rep(NA,outdim)
# 	RecTA4<-rep(NA,outdim)
# 	RecTA5<-rep(NA,outdim)
# 	RecTA6<-rep(NA,outdim)
# 	RecTA7<-rep(NA,outdim)
#   TA3.pct<-rep(NA,outdim)
#   TA4.pct<-rep(NA,outdim)
#   TA5.pct<-rep(NA,outdim)
#   TA6.pct<-rep(NA,outdim)
#   TA7.pct<-rep(NA,outdim)
# 	RecTot<-rep(NA,outdim)
# 	BYSpToSp<-rep(NA,outdim)
# 
# x1S<-X1smolt
# x2S<-X2smolt
# x3S<-X3smolt
# SmTot<-totsmolt
# OA1p<-OA1pass
# OA2p<-OA2pass
# OA3p<-OA3pass
# sadguy<-SadR
# Escape<-TotEsc
# 
# pp<-1  # a counter
# for(ss in 1:length(scen)){
#   #ssub<-scen[ss]
#   for(rr in 1:max(dat1$rep)){
#     for(g in 1:length(gmon)){ #g<-1
#       #gsub<-gmon[g]
#       for(k in 1:kmax){ #k<-1
#           for(y in minyr:(maxyr-maxage)){ #y<-8
# 
#             # get smolts by brood year and age comp
# #             BYsm1.N[pp]<-X1smolt[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+2)]
# #             BYsm2.N[pp]<-X2smolt[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+3)]
# #             BYsm3.N[pp]<-X3smolt[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+4)]
# #             BYsm1.pct[pp]<-BYsm1.N[pp]/totsmolt[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+2)]
# #             BYsm2.pct[pp]<-BYsm2.N[pp]/totsmolt[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+3)]
# #             BYsm3.pct[pp]<-BYsm3.N[pp]/totsmolt[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+4)]
# 
#               BYsm1.N[pp]<-x1S[pp+2]
#               BYsm2.N[pp]<-x2S[pp+3]
#               BYsm3.N[pp]<-x3S[pp+4]
#               BYsm1.pct[pp]<-BYsm1.N[pp]/SmTot[pp+2]
#               BYsm2.pct[pp]<-BYsm2.N[pp]/SmTot[pp+3]
#               BYsm3.pct[pp]<-BYsm3.N[pp]/SmTot[pp+4]
#             
#             # figure out what total age is for returning adults and map back to BY
#             #ocean age 1:3 fish that left as age 1 smolts
# #             BYsm1Ad1.N[pp]<-OA1pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+3)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+3)]
# #             BYsm1Ad2.N[pp]<-OA2pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+4)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+4)]
# #             BYsm1Ad3.N[pp]<-OA3pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+5)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+5)]
#               BYsm1Ad1.N[pp]<-OA1p[pp+3]*sadguy[pp+3]
#               BYsm1Ad2.N[pp]<-OA2p[pp+4]*sadguy[pp+4]
#               BYsm1Ad3.N[pp]<-OA3p[pp+5]*sadguy[pp+5]
# 
#             #ocean age 1:3 fish that left as age 2 smolts
# #             BYsm2Ad1.N[pp]<-OA1pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+4)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+4)]
# #             BYsm2Ad2.N[pp]<-OA2pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+5)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+5)]
# #             BYsm2Ad3.N[pp]<-OA3pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+6)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+6)]
#               BYsm2Ad1.N[pp]<-OA1p[pp+4]*sadguy[pp+4]
#               BYsm2Ad2.N[pp]<-OA2p[pp+5]*sadguy[pp+5]
#               BYsm2Ad3.N[pp]<-OA3p[pp+6]*sadguy[pp+6]
# 
#             #ocean age 1:3 fish that left as age 3 smolts
# #             BYsm3Ad1.N[pp]<-OA1pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+5)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+5)]
# #             BYsm3Ad2.N[pp]<-OA2pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+6)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+6)]
# #             BYsm3Ad3.N[pp]<-OA3pass[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+7)]*SadR[Gtype==gsub&Scenario==ssub&rep==rr&SubPopN==k&yr==(y+7)]
#               BYsm3Ad1.N[pp]<-OA1p[pp+5]*sadguy[pp+5]
#               BYsm3Ad2.N[pp]<-OA2p[pp+6]*sadguy[pp+6]
#               BYsm3Ad3.N[pp]<-OA3p[pp+7]*sadguy[pp+7]
# 
# 
#             #group by total age        
#             RecTA3[pp]<-BYsm1.pct[pp]*BYsm1Ad1.N[pp]
#             RecTA4[pp]<-BYsm1.pct[pp]*BYsm1Ad2.N[pp]+BYsm2.pct[pp]*BYsm2Ad1.N[pp]
#             RecTA5[pp]<-BYsm1.pct[pp]*BYsm1Ad3.N[pp]+BYsm2.pct[pp]*BYsm2Ad2.N[pp]+BYsm3.pct[pp]*BYsm3Ad1.N[pp]
#             RecTA6[pp]<-BYsm2.pct[pp]*BYsm2Ad3.N[pp]+BYsm3.pct[pp]*BYsm3Ad2.N[pp]
#             RecTA7[pp]<-BYsm3.pct[pp]*BYsm3Ad3.N[pp]
#             RecTot[pp]<-RecTA3[pp]+RecTA4[pp]+RecTA5[pp]+RecTA6[pp]+RecTA7[pp]
#             BYSpToSp[pp]<-RecTot[pp]/Escape[pp]
#             TA3.pct[pp]<-RecTA3[pp]/RecTot[pp]
#             TA4.pct[pp]<-RecTA4[pp]/RecTot[pp]
#             TA5.pct[pp]<-RecTA5[pp]/RecTot[pp]
#             TA6.pct[pp]<-RecTA6[pp]/RecTot[pp]
#             TA7.pct[pp]<-RecTA7[pp]/RecTot[pp]
# 
#             #count up
#             pp<-pp+1
#           }
#           pp<-pp+maxage
#         }
#       }
#   }
# }
# 
# dat1$OA1sp<-OA1sp
# dat1$OA2sp<-OA2sp
# dat1$OA3sp<-OA3sp
# dat1$BYsm1.N<-BYsm1.N
# dat1$BYsm2.N<-BYsm2.N
# dat1$BYsm3.N<-BYsm3.N
# dat1$BYsm1.pct<-BYsm1.pct
# dat1$BYsm2.pct<-BYsm2.pct
# dat1$BYsm3.pct<-BYsm3.pct
# dat1$BYsm1.N<-BYsm1.N
# dat1$BYsm2.N<-BYsm2.N
# dat1$BYsm3.N<-BYsm3.N
# dat1$BYsm1Ad1.N<-BYsm1Ad1.N
# dat1$BYsm1Ad2.N<-BYsm1Ad2.N
# dat1$BYsm1Ad3.N<-BYsm1Ad3.N
# dat1$BYsm2Ad1.N<-BYsm2Ad1.N
# dat1$BYsm2Ad2.N<-BYsm2Ad2.N
# dat1$BYsm2Ad3.N<-BYsm2Ad3.N
# dat1$BYsm3Ad1.N<-BYsm3Ad1.N
# dat1$BYsm3Ad2.N<-BYsm3Ad2.N
# dat1$BYsm3Ad3.N<-BYsm3Ad3.N
# dat1$RecTA3<-RecTA3
# dat1$RecTA4<-RecTA4
# dat1$RecTA5<-RecTA5
# dat1$RecTA6<-RecTA6
# dat1$RecTA7<-RecTA7
# dat1$RecTot<-RecTot
# dat1$BYSpToSp<-BYSpToSp
# dat1$TA3.pct<-TA3.pct
# dat1$TA4.pct<-TA4.pct
# dat1$TA5.pct<-TA5.pct
# dat1$TA6.pct<-TA6.pct
# dat1$TA7.pct<-TA7.pct
# 
# #setwd("C://Users//Peter//Desktop")getwd()
# write.table(dat1,
#             row.names=FALSE,sep = ",",file="MFJDWithRecruits.csv")
#############################################################################################
#############################################################################################
#############################################################################################

	# End part 1 outputs #######################################################################
	  
          
      if (which == 'all')
        {
          
          #NR=with(header, {array(rep(0,K*I*Tr*G*R),c(K,I,Tr,G,R))})
          gtype<-c("Natural","H1","N.H2","H2","N-N.H3","N-H3","N.H3-H3","N.H3-H2",
                   "H3-H2", "H3","N.H3-N.H3")
          df1col<-c("MCrun","year","site.num","site.name","stage","g.num","g.name","Ni")
          df1<-array(NA,c(header$K*header$I*header$Tr*header$G*header$R+1,8))
          df1[1,]<-df1col
          ii<-2
          for(a in 1:header$R )
          {
            for(b in 1:header$K)
            {
              for(c in 1:header$Tr)
              {
                for(d in 1:header$G)
                {
                  for(e in 1:header$I)
                  {
                    df1[ii,]<-c(a,c,b,header$Site.Names[b],e,d,gtype[d],NR[b,e,c,d,a])
                    ii<-ii+1
                  }
                }
              }
            }
          }
          #names(df1)<-df1col
          #write.csv(df1,"NI.csv",col.names=NA,row.names=NA)
          Fname<-paste("jdout_NR_R",Sys.time(),".csv")
          Fname<-gsub("-","",Fname)
          Fname<-gsub(" ","_",Fname)
          Fname<-gsub(":","",Fname)
          write.table(df1, file = Fname,row.names=FALSE, na="",col.names=FALSE, sep=",")
          
          
          
          #NR=with(header, {array(rep(0,K*I*Tr*G*R),c(K,I,Tr,G,R))})
          gtype<-c("Natural","H1","N.H2","H2","N-N.H3","N-H3","N.H3-H3","N.H3-H2",
                   "H3-H2", "H3","N.H3-N.H3")
          df2col<-c("MCrun","year","site.num","site.name","stage","g.num","g.name","Ni5")
          df2<-array(NA,c(header$K*header$I*header$Tr*header$G*header$R+1,8))
          df2[1,]<-df2col
          ii<-2
          for(a in 1:header$R )
          {
            for(b in 1:header$K)
            {
              for(c in 1:header$Tr)
              {
                for(d in 1:header$G)
                {
                  for(e in 1:header$I5)
                  {
                    df2[ii,]<-c(a,c,b,header$Site.Names[b],e,d,gtype[d],N5R[b,e,c,d,a])
                    ii<-ii+1
                  }
                }
              }
            }
          }
          #names(df1)<-df1col
          #write.csv(df1,"NI.csv",col.names=NA,row.names=NA)
          Fname<-paste("jdout_N5R_R",Sys.time(),".csv")
          Fname<-gsub("-","",Fname)
          Fname<-gsub(" ","_",Fname)
          Fname<-gsub(":","",Fname)
          write.table(df2, file = Fname,row.names=FALSE, na="",col.names=FALSE, sep=",")
          
          
          #NR=with(header, {array(rep(0,K*I*Tr*G*R),c(K,I,Tr,G,R))})
          gtype<-c("Natural","H1","N.H2","H2","N-N.H3","N-H3","N.H3-H3","N.H3-H2",
                   "H3-H2", "H3","N.H3-N.H3")
          df3col<-c("MCrun","year","site.num","site.name","g.num","g.name","Spawners")
          df3<-array(NA,c(header$K*header$Tr*header$G*header$R+1,7))
          df3[1,]<-df3col
          ii<-2
          for(a in 1:header$R )
          {
            for(b in 1:header$K)
            {
              for(c in 1:header$Tr)
              {
                for(d in 1:header$G)
                {
                  df3[ii,]<-c(a,c,b,header$Site.Names[b],d,gtype[d],Spawners_NR[b,c,d,a])
                  ii<-ii+1
                }
              }
            }
          }
          #names(df1)<-df1col
          #write.csv(df1,"NI.csv",col.names=NA,row.names=NA)
          Fname<-paste("jdout_Spawners_NR_R",Sys.time(),".csv")
          Fname<-gsub("-","",Fname)
          Fname<-gsub(" ","_",Fname)
          Fname<-gsub(":","",Fname)
          write.table(df3, file = Fname,row.names=FALSE, na="",col.names=FALSE, sep=",") 
          
          
        }

      
      
      
      
      
      
  }
