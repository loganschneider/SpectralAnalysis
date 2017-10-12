#This function requires the assessment period folder name as input
#There must be a "period"_analysis folder from running the NSRR SpectralTrainFig MATLAB script
#It then automates the sleep feature definition/extraction (ultradian cycles, lights off, etc.)
spectranal <- function(period="Baseline") {
    setwd("D:/")
    edflist <- list.files(paste(getwd(),"/",period,"/",sep=""),pattern = ".edf$")
    ###Had to get rid of the messy parts of the file names (e.g. "-12(BL))
    #list <- gsub(pattern = "\\.edf$", "", edflist)
    list <- substr(edflist,start=1,stop=7)
    
    library("xlsx")
    library("data.table")
    
    speclist <- list.files(paste(period,"_analysis/",sep=""),pattern = ".spectral.xlsx$")
    setwd(paste(period,"_analysis/",sep=""))
    for(i in list) {
        specanal <- read.xlsx(grep(pattern = paste("^",i,".*detail.spectral.xlsx$",sep=""),speclist,value = T),
                              sheetName = "Sheet1",
                              startRow = 2,
                              header = T)
        names(specanal) <- c("epoch",
                             "charHypno",
                             "stage",
                             "PD0",
                             "PD0.5",
                             "PD1",
                             "PD1.5",
                             "PD2",
                             "PD2.5",
                             "PD3",
                             "PD3.5",
                             "PD4",
                             "PD4.5",
                             "PD5",
                             "PD5.5",
                             "PD6",
                             "PD6.5",
                             "PD7",
                             "PD7.5",
                             "PD8",
                             "PD8.5",
                             "PD9",
                             "PD9.5",
                             "PD10",
                             "PD10.5",
                             "PD11",
                             "PD11.5",
                             "PD12",
                             "PD12.5",
                             "PD13",
                             "PD13.5",
                             "PD14",
                             "PD14.5",
                             "PD15",
                             "PD15.5",
                             "PD16",
                             "PD16.5",
                             "PD17",
                             "PD17.5",
                             "PD18",
                             "PD18.5",
                             "PD19",
                             "PD19.5",
                             "PD20",
                             "PD20.5",
                             "PD21",
                             "PD21.5",
                             "PD22",
                             "PD22.5",
                             "PD23",
                             "PD23.5",
                             "PD24",
                             "PD24.5",
                             "PD25",
                             "P0",
                             "P0.5",
                             "P1",
                             "P1.5",
                             "P2",
                             "P2.5",
                             "P3",
                             "P3.5",
                             "P4",
                             "P4.5",
                             "P5",
                             "P5.5",
                             "P6",
                             "P6.5",
                             "P7",
                             "P7.5",
                             "P8",
                             "P8.5",
                             "P9",
                             "P9.5",
                             "P10",
                             "P10.5",
                             "P11",
                             "P11.5",
                             "P12",
                             "P12.5",
                             "P13",
                             "P13.5",
                             "P14",
                             "P14.5",
                             "P15",
                             "P15.5",
                             "P16",
                             "P16.5",
                             "P17",
                             "P17.5",
                             "P18",
                             "P18.5",
                             "P19",
                             "P19.5",
                             "P20",
                             "P20.5",
                             "P21",
                             "P21.5",
                             "P22",
                             "P22.5",
                             "P23",
                             "P23.5",
                             "P24",
                             "P24.5",
                             "P25",
                             "sleepMask",
                             "nremMask",
                             "remMask",
                             "deltaArtifactMask",
                             "betaArtifactMask",
                             "artifactMask")
        firststage <- specanal[match(unique(specanal$stage), specanal$stage),c("epoch","stage")]
        LOff <- firststage[which(firststage$stage == 0),"epoch"] #Lights off
        SO <- min(firststage[which(firststage$stage == 1),"epoch"],
                  firststage[which(firststage$stage == 2),"epoch"],
                  firststage[which(firststage$stage == 3),"epoch"],
                  firststage[which(firststage$stage == 4),"epoch"],
                  firststage[which(firststage$stage == 5),"epoch"]) #Sleep onset by first of any stage
        SOL <- round((SO-LOff)/2,1) #Sleep onset latency in min
        RL <- round((firststage[which(firststage$stage == 5),"epoch"]-SO)/2,1) #REM latency in min
        Final <- max(specanal[which(specanal$stage == 1),"epoch"],
                   specanal[which(specanal$stage == 2),"epoch"],
                   specanal[which(specanal$stage == 3),"epoch"],
                   specanal[which(specanal$stage == 4),"epoch"],
                   specanal[which(specanal$stage == 5),"epoch"]) #Final sleep epoch
        LOn <- max(specanal[which(specanal$stage == 0),"epoch"]) #Lights on...may technically be erroneous (not sure of the true marker)
        
        #Create a data frame for epochs of interest
        markers <- data.frame("Marker"=c("Lights off","Sleep onset","Final sleep epoch","Lights on"),
                              "epoch"=c(LOff, SO, Final, LOn))
        
        #Now working on stages and cycles
        specsub <- subset.data.frame(specanal,specanal$stage<6)
        id <- cumsum(c(TRUE,diff(specsub$stage)!=0)) #number of stage transitions
        #Identify each stage bout (even single-epoch instances) for ultradian calculation
        stagestats <- as.data.frame(cbind(stage=tapply(specsub$stage,id,max),
              firstEpoch=tapply(specsub$epoch,id,min),
              lastEpoch=tapply(specsub$epoch,id,max)))
        #Calculation stage time from epochs
        for(j in 1:nrow(stagestats)) {
            stagestats[j,"DurationInMin"] <- round((stagestats[j,"lastEpoch"]-stagestats[j,"firstEpoch"]+1)/2,1)
        }
        
        #Define each bout of REM (to be able to determine cycles by ultradian bouts of REM)
        stagestats$REMblock <- ifelse(stagestats$stage==5,1,0)
        dt <- as.data.table(stagestats)
        
        #Sequentially number each bout of REMS
        dt[, REMbout := 1:(.N), by="REMblock"]
        stagestats <- as.data.frame(dt)
        stagestats$REMbout[stagestats$REMblock == 0] <- 0
        
        #Define number of cycles of REMS for demarcation of ultradian cycles
        REMs <- stagestats[which(stagestats$REMblock==1),]; REMs #Just to check
        cycle=1; REMs[1,"uCycle"] <- cycle
        cycledf <- data.frame("CycleNum"=as.numeric(1), "StartEpoch"=as.numeric(SO), "EndEpoch"=numeric(1))
        for(k in 2:nrow(REMs)) {
            #This conditional starts things off in the event of only 1 REM cycle
            if(k==nrow(REMs) & is.na(REMs[nrow(REMs),"uCycle"])) { 
                cycledf[cycle,"EndEpoch"] <- REMs[k,"lastEpoch"]
            } else if(nrow(REMs)==1 && length(RL)!=0) {
                #Needed this addition because of cases like P100407 with only 1 bout of REM
                cycledf[cycle,"EndEpoch"] <- REMs[k-1,"lastEpoch"]
                break #Loop is essentially done, as only 1 REM bout
            } else if(nrow(REMs)==1 && length(RL)==0){
                #Defining single ultradian cycle as whole sleep period if no REM achieved
                cycledf[cycle,"EndEpoch"] <- Final
                break #If no REM achieved, no need to loop
            }
            if(REMs[k,"firstEpoch"]<(REMs[k-1,"lastEpoch"]+90)) {
                #If next REM bout is <45 min away, then include it in this REM cycle
                REMs[k,"uCycle"] <- cycle
            } else {
                #If next REM bout is >45 min away, then start a new ultradian/REM cycle
                cycledf[cycle,"EndEpoch"] <- REMs[k-1,"lastEpoch"]
                cycle=cycle+1
                REMs[k,"uCycle"] <- cycle
                cycledf[cycle,"CycleNum"] <- cycle
                cycledf[cycle,"StartEpoch"] <- REMs[k-1,"lastEpoch"]+1
            }
            #This conditional wraps things up by designating next ultradian/REM cycle
            if(k==nrow(REMs) & is.na(REMs[nrow(REMs),"uCycle"])) { 
                REMs[k,"uCycle"] <- cycle
                cycledf[cycle,"EndEpoch"] <- REMs[k,"lastEpoch"]
            } else if(k==nrow(REMs)) { 
                REMs[k,"uCycle"] <- cycle
                cycledf[cycle,"EndEpoch"] <- REMs[k,"lastEpoch"]
            }
        }
        
        #Generate a total list from first sleep epoch to last sleep epoch of cycles
        #This will be used for creating an ultradian cycle marker in specanal
        cycling <- data.frame("StartStop"=numeric(0))
        for(l in 1:nrow(cycledf)) {
            cycling[2*l-1,"StartStop"] <- cycledf[l,"StartEpoch"]
            cycling[2*l,"StartStop"] <- cycledf[l,"EndEpoch"]
        }
        
        #Necessary to add a final cycle and epoch, if not having terminal awakening from REM
        if(cycling[nrow(cycling),"StartStop"] != Final) {
            cycling[nrow(cycling)+1,"StartStop"] <- cycling[nrow(cycling),"StartStop"]+1
            cycling[nrow(cycling)+1,"StartStop"] <- Final
        }
        
        #Add ultradian cycle ("uCycle") information to specanal
        for(m in 1:(nrow(cycling)/2)) {
            specanal[cycling[2*m-1,"StartStop"]:cycling[2*m,"StartStop"],"uCycle"] <- m
        }
        
        #Write out the tables for analysis
        write.table(markers,file=paste(i,"markers.txt",sep="_"),quote=F,row.names=F)
        write.table(stagestats,file=paste(i,"stagestats.txt",sep="_"),quote=F,row.names=F)
        write.table(REMs,file=paste(i,"REMbouts.txt",sep="_"),quote=F,row.names=F)
        write.table(cycledf,file=paste(i,"ultradiancycles.txt",sep="_"),quote=F,row.names=F)
        write.csv(specanal,file=paste(i,".spectral.csv",sep=""))
        
        success <- paste("Successfully completed",i,sep=" ")
        write.table(success,"Last successful subject.txt",row.names = F, col.names = F)
    }
}
