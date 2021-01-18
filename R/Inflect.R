#' This function analyzes raw abundance data from a Thermal Profiling experiment and calculates melt temperatures and melt shifts for each protein in the experiment.
#' An example line of code that could be run to test the program function is shown below.
#' The directory below is an example directory but should be replaced with the path that the source files are present.
#' This directory is also the location where the output files including pdf and Excel will be transferred once the calculations are complete.
#' The source files should be Excel format or xlsx and there should be two files for each replicate condition.
#' For example if there is one replicate, there should be one file titled Condition 1.xlsx and another file titled Control 1.xlsx. Data from the Control or vehicle will be in the Control file and data from the condition or treatment will be in the Condition file.
#' Each spreadsheet file should have the format where the first column header should be Accession and all of the values in this column are accession numbers for the proteins.
#' There should be a column then for each temperature at which the experiment was run. For example if there are 8 temperatures, there should be 8 columns after the accession column.
#' Each temperature column can have its own label but the values in each cell should correspond to the abundance values for each protein and at each temperature.
#' directory<-"/Users/folder_1/folder_2"
#' Temperature<-c(35,45,50,55,60,75)
#' Rsq<-0.95
#' NumSD<-2
#' NReps<-3
#' Inflect(directory,Temperature,Rsq,NumSD,NReps)
#' @param directory the directory where the source data files to be analyzed are saved. This is also the location where the results will be saved.
#' @param Temperature the temperatures from the heat treatment procedure. An example entry Temperature<-c(25,35,39.3,50.1,55.2,60.7,74.9,90)
#' @param Rsq the cutoff to be used for the melt shift curve fit. An example entry would be 0.95
#' @param NumSD the standard deviation cutoff to be used for the calculated melt shifts. For example, if NumSD = 2, proteins with melt shifts greater than 2 standard deviations from the mean will be considered significant.
#' @param NReps the number of replicate experiments to be analyzed
#' @importFrom readxl read_excel
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
#' @return xlsx files with calculated melt shift for each protein in the experiment
Inflect<-function(directory,Temperature,Rsq,NumSD,NReps){

SourcePath<-directory

#This section of code repeats the analysis function for every replicate data set.
RepNum<-1
repeat{
  MasterDirectory<-paste(directory,paste("Rep",RepNum),sep="/")
  dir.create(MasterDirectory)
  dir.create(paste(MasterDirectory,"Curves",sep="/"))
  dir.create(paste(MasterDirectory,"SigCurves",sep="/"))
  OutputPath_Curves=paste(MasterDirectory,"Curves",sep="/")
  OutputPath_SigCurves=paste(MasterDirectory,"SigCurves",sep="/")
  OutputPath =MasterDirectory
  Inflect<-InflectWorkflow(Rsq,NumSD,Temperature,RepNum,SourcePath,OutputPath)
  RepNum<-RepNum+1
  if (RepNum > NReps){
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
}

