#' This function analyzes raw abundance data from a Thermal Proteome Profiling experiment and calculates melt temperatures and melt shifts for each protein in the experiment.
#' @param directory the directory where the source data files to be analyzed are saved. This is also the location where the results will be saved.
#' @param Temperature the temperatures from the heat treatment procedure. An example entry Temperature<-c(25,35,39.3,50.1,55.2,60.7,74.9,90)
#' @param Rsq the cutoff to be used for the melt shift curve fit. An example entry would be 0.95
#' @param NumSD the standard deviation cutoff to be used for the calculated melt shifts. For example, if NumSD = 2, proteins with melt shifts greater than 2 standard deviations from the mean will be considered significant.
#' @param NReps the number of replicate experiments to be analyzed
#' @importFrom readxl read_excel
#' @importFrom readxl read_xlsx
#' @importFrom writexl write_xlsx
#' @importFrom stats median
#' @importFrom stats nls
#' @importFrom stats sd
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom tidyr spread
#' @importFrom grDevices dev.off
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics axis
#' @importFrom graphics legend
#' @importFrom plotrix addtable2plot
#' @importFrom data.table data.table
#' @importFrom grDevices pdf
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_set
#' @importFrom ggplot2 theme_gray
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 ggsave
#' @importFrom UpSetR upset
#' @importFrom UpSetR fromList
#' @return xlsx files with calculated melt shift for each protein in the experiment along with Upset plots that show the overlap in number of proteins stabilized and destabilized between each replicate
Inflect<-function(directory,Temperature,Rsq,NumSD,NReps){

SourcePath<-directory

#This section of code repeats the analysis function for every replicate data set.
Rep<-1
repeat{
  MasterDirectory<-paste(directory,paste("Rep",Rep),sep="/")
  dir.create(MasterDirectory)
  dir.create(paste(MasterDirectory,"Curves",sep="/"))
  dir.create(paste(MasterDirectory,"SigCurves",sep="/"))
  OutputPath_Curves=paste(MasterDirectory,"Curves",sep="/")
  OutputPath_SigCurves=paste(MasterDirectory,"SigCurves",sep="/")
  OutputPath =MasterDirectory
  Inflect<-InflectWorkflow(Rsq,NumSD,Temperature,Rep,SourcePath,OutputPath)
  Rep<-Rep+1
  if (Rep > NReps){
    break
  }
}

#This section of code summarizes the results from the replicates by combining all of the significant proteins and their melt shifts into a single table. This section also separates destabilized and stabilized proteins so that they can be analyzed by Upset plot.
Rep<-1
DataRows<-0
repeat{
  DataRows<-nrow(read_excel(file.path(paste(paste(SourcePath,paste("Rep",Rep),sep="/")),"SignificantResults.xlsx")))+DataRows
  Rep<-Rep+1
  if (Rep > NReps){
    break
  }
}

MeltShiftDataSubset = matrix(NA,nrow=DataRows ,ncol=3)
MeltShiftDataUpset = as.data.frame(matrix(NA,nrow=DataRows ,ncol=NReps))
MeltShiftDataUpsetDS = as.data.frame(matrix(NA,nrow=DataRows ,ncol=NReps*2))
colnames(MeltShiftDataSubset)<-c("Replicate","Accession","MeltShift")
Rep<-1
Row<-0
repeat{
  OutputPath =directory
  MeltShiftData<-read_excel(file.path(paste(paste(SourcePath,paste("Rep",Rep),sep="/")),"SignificantResults.xlsx"))
  Rowcount<-as.numeric(nrow(MeltShiftData))
  MeltShiftDataSubset[(Row+1):(Row+Rowcount),1]<-paste("Replicate",Rep)
  MeltShiftDataSubset[(Row+1):(Row+Rowcount),2]<-MeltShiftData$Accession
  MeltShiftDataSubset[(Row+1):(Row+Rowcount),3]<-as.numeric(MeltShiftData$DeltaTmControlMinusCondition)
  MeltShiftDataUpset[1:Rowcount,Rep]<-MeltShiftData$Accession
  colnames(MeltShiftDataUpset)[Rep]<-paste("Rep",Rep)
  RowcountDestab<-as.numeric(nrow(MeltShiftData[MeltShiftData$DeltaTmControlMinusCondition>0,]))
  RowcountStab<-as.numeric(nrow(MeltShiftData[MeltShiftData$DeltaTmControlMinusCondition<0,]))
  MeltShiftDataUpsetDS[1:RowcountDestab,Rep]<-MeltShiftData[MeltShiftData$DeltaTmControlMinusCondition>0,2]
  colnames(MeltShiftDataUpsetDS)[Rep]<-paste("Destablized Rep",Rep)
  MeltShiftDataUpsetDS[1:RowcountStab,Rep+NReps]<-MeltShiftData[MeltShiftData$DeltaTmControlMinusCondition<0,2]
  colnames(MeltShiftDataUpsetDS)[Rep+NReps]<-paste("Stabilized Rep",Rep)
  Row<-Row+Rowcount
  Rep<-Rep+1
  if (Rep > NReps){
    break
  }
}
MeltShiftDataSubsetWide <- spread(as.data.frame(MeltShiftDataSubset), "Replicate", "MeltShift")
MeltShiftDataSubsetStabilized = matrix(NA,nrow=DataRows ,ncol=3)
MeltShiftDataSubsetDestabilized = matrix(NA,nrow=DataRows ,ncol=3)
colnames(MeltShiftDataSubsetStabilized)<-c("Replicate","Accession","MeltShift")
colnames(MeltShiftDataSubsetDestabilized)<-c("Replicate","Accession","MeltShift")
MeltShiftDataSubset<-as.data.frame(MeltShiftDataSubset)
MeltShiftDataSubsetStabilized<-MeltShiftDataSubset[MeltShiftDataSubset$MeltShift<0,]
MeltShiftDataSubsetDestabilized<-MeltShiftDataSubset[MeltShiftDataSubset$MeltShift>0,]
MeltShiftDataSubsetStabilizedWide <- spread(as.data.frame(MeltShiftDataSubsetStabilized), "Replicate", "MeltShift")
MeltShiftDataSubsetDestabilizedWide <- spread(as.data.frame(MeltShiftDataSubsetDestabilized), "Replicate", "MeltShift")

#This section of code summarizes results in a table that is added to the directory along with files that can be used for analysis by other programs.
write_xlsx(MeltShiftDataSubsetWide,paste(OutputPath,"SummaryResults.xlsx",sep="/"))
write_xlsx(MeltShiftDataUpset,paste(OutputPath,"AllSignificant.xlsx",sep="/"))
write_xlsx(MeltShiftDataUpsetDS,paste(OutputPath,"AllSignificantStabDestab.xlsx",sep="/"))
DestabList<-MeltShiftDataUpsetDS[c(1:NReps)]
StabList<-MeltShiftDataUpsetDS[c((NReps+1):(NReps*2))]

write_xlsx(DestabList,paste(OutputPath,"AllSignificantDestabilized.xlsx",sep="/"))
write_xlsx(StabList,paste(OutputPath,"AllSignificantStabilized.xlsx",sep="/"))

DestabFile<-as.data.frame(read_excel(paste(directory,"AllSignificantDestabilized.xlsx",sep="/")))
StabFile<-as.data.frame(read_excel(paste(directory,"AllSignificantStabilized.xlsx",sep="/")))
pdf(file=paste(directory,"DestabilizedUpset.pdf",sep="/"),onefile = FALSE)
print(upset(fromList(DestabFile),nsets=NReps))
dev.off()
pdf(file=paste(directory,"StabilizedUpset.pdf",sep="/"),onefile = FALSE)
print(upset(fromList(StabFile),nsets=NReps))
dev.off()


}

