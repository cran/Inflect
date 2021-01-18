#' This function analyzes raw abundance data from a Thermal Profiling experiment and calculates melt temperatures and melt shifts for each protein in the experiment.
#' See the Inflect function file for an example of program that could be executed to test the program function.
#' @param Temperature the temperatures from the heat treatment procedure. An example entry Temperature<-c(25,35,39.3,50.1,55.2,60.7,74.9,90)
#' @param Rsq the cutoff to be used for the melt shift curve fit. An example entry would be 0.95
#' @param NumSD the standard deviation cutoff to be used for the calculated melt shifts. For example, if NumSD = 2, proteins with melt shifts greater than 2 standard deviations from the mean will be considered significant.
#' @param Rep the replicate number that is being analyzed
#' @param SourcePath The path for the source data
#' @param OutputPath The path for the output data
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
InflectWorkflow<-function(Rsq,NumSD,Temperature,Rep,SourcePath,OutputPath){

  OutputPath_Curves=paste(OutputPath,"Curves",sep="/")
  OutputPath_SigCurves=paste(OutputPath,"SigCurves",sep="/")

  #This chunk prepares data for log fitting.

  FileCondition<-paste(paste("Condition",Rep),"xlsx",sep=".")
  FileControl<-paste(paste("Control",Rep),"xlsx",sep=".")

  NumberTemperatures<-as.numeric(NROW(Temperature))
  ConditionData <- as.data.frame(read_excel(file.path(paste(SourcePath,FileCondition,sep="/"))))
  ControlData <- as.data.frame(read_excel(file.path(paste(SourcePath,FileControl,sep="/"))))
  Data_Control <- ControlData[,c(2:(1+NumberTemperatures))]
  Data_Condition <- ConditionData[,c(2:(1+NumberTemperatures))]
  Protein_Control <- ControlData[c(1)]
  Protein_Condition <- ConditionData[c(1)]
  Data_Control_Norm<-Data_Control[,1:ncol(Data_Control)]/Data_Control[,1]
  Data_Condition_Norm<- Data_Condition[,1:ncol(Data_Condition)]/Data_Condition[,1]
  All_Control_Norm<-data.frame(Protein_Control,Data_Control_Norm)
  All_Condition_Norm<-data.frame(Protein_Condition,Data_Condition_Norm)

  #Omits rows with N/A values and writes to file the proteins at this point
  All_Control_Norm_Omit<-na.omit(All_Control_Norm)
  All_Condition_Norm_Omit<-na.omit(All_Condition_Norm)
  Proteins_Control_Norm_Omit<-All_Control_Norm_Omit[1]
  Proteins_Condition_Norm_Omit<-All_Condition_Norm_Omit[1]
  Data_Control_Norm_Omit<-((All_Control_Norm_Omit[2:(1+NumberTemperatures)]))
  Data_Condition_Norm_Omit<-(All_Condition_Norm_Omit[2:(1+NumberTemperatures)])
  write_xlsx(Data_Control_Norm_Omit,paste(OutputPath,"NormalizedControlResults.xlsx",sep="/"))
  write_xlsx(Data_Condition_Norm_Omit,paste(OutputPath,"NormalizedConditionResults.xlsx",sep="/"))

  #Calculates the median fold change abundance from each temperature
  Data_Control_Norm_Omit_Median<-apply(Data_Control_Norm_Omit, 2, FUN=median)
  Data_Condition_Norm_Omit_Median<-apply(Data_Condition_Norm_Omit, 2, FUN=median)

  ControlMedian<-Data_Control_Norm_Omit_Median[1:NumberTemperatures]
  ConditionMedian<-Data_Condition_Norm_Omit_Median[1:NumberTemperatures]

  Proteins_Control_Norm_Omit<-All_Control_Norm_Omit[c(1)]
  Proteins_Condition_Norm_Omit<-All_Condition_Norm_Omit[c(1)]

  #This portion of code calculates the 4 parameter logistic curve for the normalized values at each temperature.

  ControlNormBothCorrect<-FPLFit_Correction(ControlMedian,Data_Control_Norm_Omit,"Control",Temperature)
  ConditionNormBothCorrect<-FPLFit_Correction(ConditionMedian,Data_Condition_Norm_Omit,"Condition",Temperature)

  #This section of code calculates the 4 parameter logistic variables for each of the proteins from both Control and Condition conditions.
  Data_Norm_Omit<-Data_Control_Norm_Omit
  NormBothCorrect<-ControlNormBothCorrect

  DataParametersControl<-FPLFit(Data_Control_Norm_Omit,ControlNormBothCorrect,"Control",Temperature,NumberTemperatures)
  DataParametersCondition<-FPLFit(Data_Condition_Norm_Omit,ConditionNormBothCorrect,"Condition",Temperature,NumberTemperatures)

  #This section of code isolates proteins that have a Rsquare of a certain value

  DataParametersControlResults <- data.frame(Proteins_Control_Norm_Omit, DataParametersControl, ControlNormBothCorrect)
  DataParametersConditionResults <- data.frame(Proteins_Condition_Norm_Omit, DataParametersCondition, ConditionNormBothCorrect)
  DataParametersResultsAll <-merge(DataParametersControlResults,DataParametersConditionResults,by="Accession")

  significance_condition <- DataParametersResultsAll[c(6)] >= Rsq & DataParametersResultsAll[c(NumberTemperatures+11)] >= Rsq

  SignificantAll<- DataParametersResultsAll[significance_condition, c(1:NCOL(DataParametersResultsAll))]
  ProteinsFilter<- DataParametersResultsAll[significance_condition, c(1)]

  DeltaTmControlMinusCondition<-SignificantAll[,3]-SignificantAll[,8+NumberTemperatures]
  DeltaTmDataSet<-data.table(SignificantAll,DeltaTmControlMinusCondition)
  DeltaTmSort<-DeltaTmDataSet[order(-DeltaTmControlMinusCondition)]
  Observation <- 1:NROW(ProteinsFilter)
  DataParametersAnalysisResults <- data.frame(Observation,DeltaTmSort)

  Observation <- as.numeric(DataParametersAnalysisResults[,1])
  DeltTm <- DataParametersAnalysisResults[,NCOL(DataParametersAnalysisResults)]
  TmData <- data.table(Observation, DeltTm)

  SDH<-mean(DeltTm)+NumSD*sd(DeltTm)
  SDL<-mean(DeltTm)-NumSD*sd(DeltTm)

  significance_condition <- DataParametersAnalysisResults[c(NCOL(DataParametersAnalysisResults))] >= SDH | DataParametersAnalysisResults[c(NCOL(DataParametersAnalysisResults))] <= SDL
  SignificantAll_SD<- DataParametersAnalysisResults[significance_condition, c(1:NCOL(DataParametersAnalysisResults))]

  #This section of code creates the output files from the analysis including pdf files each with a prefix corresponding to the accession number or protein name in the original source file.

  write_xlsx(DataParametersResultsAll,paste(OutputPath,"Results.xlsx",sep="/"))
  write_xlsx(SignificantAll_SD,paste(OutputPath,"SignificantResults.xlsx",sep="/"))

  n <- 1
  repeat{
    pdf(paste(OutputPath_Curves,paste(DataParametersAnalysisResults[n,2],"pdf",sep="."),sep="/"))
    plot(x = Temperature, y = DataParametersAnalysisResults[n,c(8:(NumberTemperatures+7))],pch = 2, frame = TRUE,xlab = "Temperature (C)", ylab = "Normalized Abundance",col = "blue",axes=FALSE,ylim=c(0,1.5),xlim=c(20,100),main=DataParametersAnalysisResults[n,2])
    Temperaturevals<-seq(min(Temperature), max(Temperature), by=1)

    lines(Temperaturevals ,(DataParametersAnalysisResults[n,c(5)]+((DataParametersAnalysisResults[n,c(6)]-DataParametersAnalysisResults[n,c(5)])/(1+exp(-DataParametersAnalysisResults[n,c(3)]*(Temperaturevals-DataParametersAnalysisResults[n,c(4)]))))),col="blue")

    points(Temperature, DataParametersAnalysisResults[n,c((13+NumberTemperatures):(12+(2*NumberTemperatures)))], col="red", pch=8,lty=1)

    lines(Temperaturevals<-seq(min(Temperature), max(Temperature), by=1),(DataParametersAnalysisResults[n,c(NumberTemperatures+10)]+((DataParametersAnalysisResults[n,c(NumberTemperatures+11)]-DataParametersAnalysisResults[n,c(NumberTemperatures+10)])/(1+exp(-DataParametersAnalysisResults[n,c(NumberTemperatures+8)]*(Temperaturevals-DataParametersAnalysisResults[n,c(NumberTemperatures+9)]))))),lty=2,col="red")

    axis(side=1, at=seq(20,100,by=5))
    axis(side=2, at=seq(0,1.5, by=.1))
    legend(60,1.3,legend=c("Control","Condition"), col=c("blue","red"),pch=c(2,8),lty=c(1,2))
    plottable=matrix(data=NA,nrow=3,ncol=2)
    colnames(plottable)<-c("Condition","Temperature")
    plottable[1,1]<-"Control Tm"
    plottable[2,1]<-"Condition Tm"
    plottable[3,1]<-"Delta Tm"
    plottable[1,2]<-round(DataParametersAnalysisResults[n,c(4)],2)
    plottable[2,2]<-round(DataParametersAnalysisResults[n,c(NumberTemperatures+9)],2)
    plottable[3,2]<-round((DataParametersAnalysisResults[n,c(4)])-(DataParametersAnalysisResults[n,c(NumberTemperatures+9)]),2)
    addtable2plot(60,0.8,plottable,hlines=TRUE,vlines=TRUE,bty="o",bg="gray")
    dev.off()
    n=n+1
    if (n > NROW(DataParametersAnalysisResults)){
      break
    }
  }

  n <- 1
  repeat{
    pdf(paste(OutputPath_SigCurves,paste(SignificantAll_SD[n,2],"pdf",sep="."),sep="/"))
    plot(x = Temperature, y = SignificantAll_SD[n,c(8:(NumberTemperatures+7))],pch = 2, frame = TRUE,xlab = "Temperature (C)", ylab = "Normalized Abundance",col = "blue",axes=FALSE,ylim=c(0,1.5),xlim=c(20,100),main=SignificantAll_SD[n,2])
    Temperaturevals<-seq(min(Temperature), max(Temperature), by=1)

    lines(Temperaturevals ,(SignificantAll_SD[n,c(5)]+((SignificantAll_SD[n,c(6)]-SignificantAll_SD[n,c(5)])/(1+exp(-SignificantAll_SD[n,c(3)]*(Temperaturevals-SignificantAll_SD[n,c(4)]))))),col="blue")

    points(Temperature, SignificantAll_SD[n,c((13+NumberTemperatures):(12+(2*NumberTemperatures)))], col="red", pch=8,lty=1)

    lines(Temperaturevals<-seq(min(Temperature), max(Temperature), by=1),(SignificantAll_SD[n,c(NumberTemperatures+10)]+((SignificantAll_SD[n,c(NumberTemperatures+11)]-SignificantAll_SD[n,c(NumberTemperatures+10)])/(1+exp(-SignificantAll_SD[n,c(NumberTemperatures+8)]*(Temperaturevals-SignificantAll_SD[n,c(NumberTemperatures+9)]))))),lty=2,col="red")

    axis(side=1, at=seq(20,100,by=5))
    axis(side=2, at=seq(0,1.5, by=.1))
    legend(65,1.3,legend=c("Control","Condition"), col=c("blue","red"),pch=c(2,8),lty=c(1,2))
    plottable=matrix(data=NA,nrow=3,ncol=2)
    colnames(plottable)<-c("Condition","Temperature")
    plottable[1,1]<-"Control Tm"
    plottable[2,1]<-"Condition Tm"
    plottable[3,1]<-"Delta Tm"
    plottable[1,2]<-round(SignificantAll_SD[n,c(4)],2)
    plottable[2,2]<-round(SignificantAll_SD[n,c(NumberTemperatures+9)],2)
    plottable[3,2]<-round((SignificantAll_SD[n,c(4)])-(SignificantAll_SD[n,c(NumberTemperatures+9)]),2)
    addtable2plot(65,0.8,plottable,hlines=TRUE,vlines=TRUE,bty="o",bg="gray")

    dev.off()
    n=n+1
    if (n > NROW(SignificantAll_SD)){
      break
    }
  }

  WaterfallPlot <- ggplot(DataParametersAnalysisResults[c(1,NCOL(DataParametersAnalysisResults))],
                          aes(x = Observation,y = DeltTm)) +
    geom_point(size=3) +
    labs(x = "Rank of Melting Point Difference Control and Condition",
         y = "Delta Tm (C)") +
    ylim(c(1.1*min(DataParametersAnalysisResults$DeltaTmControlMinusCondition),1.1*max(DataParametersAnalysisResults$DeltaTmControlMinusCondition))) +
    theme_minimal() +
    theme(legend.position="none") +theme_set(theme_gray(base_size = 15))+
    geom_hline(yintercept = SDH, linetype = "dashed", color = "orange") +
    geom_hline(yintercept=SDL, linetype = "dashed", color = "orange")

  plot(WaterfallPlot)

  ggsave(filename = "WaterfallPlot.pdf",
         path = OutputPath_Curves,
         plot = WaterfallPlot,
         device = "pdf",
         width = 10,
         height = 10,
         units = "in",
         useDingbats = FALSE)

}
