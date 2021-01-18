#' This function determines the 4 parameter log fit constants for all proteins in the control and condition files.
#' See the Inflect function file for an example of program that could be executed to test the program function.
#' @param Median Median abundance values
#' @param Data_Norm_Omit This is the normalized data where data that does not meet criteria has been omitted.
#' @param Condition This is the control or condition that is being analyzed
#' @param Temperature The temperatures from the heat treatment
#' @return set of correction constants for each control and condition
FPLFit_Correction<-function(Median,Data_Norm_Omit,Condition,Temperature){
  TotalCellsPer<-as.numeric(NROW(Data_Norm_Omit)*NCOL(Data_Norm_Omit))
  TotalRowsPer<-as.numeric(NROW(Data_Norm_Omit))
  TotalColsPer<-as.numeric(NCOL(Data_Norm_Omit))
  df <- data.frame(Temperature, Median)
  CurveFit <- nls(Median ~ I(c+((d-c)/(1+exp(-a*(Temperature-b))))), data = df, start=list(a=-1,b=mean(Temperature),c=min(Median),d=max(Median)), trace = F)
  CCurvea1<-summary(CurveFit)
  CCurvea2<-CCurvea1$coefficients
  CCurvea<-CCurvea2[1,1]
  CCurveb<-CCurvea2[2,1]
  CCurvec<-CCurvea2[3,1]
  CCurved<-CCurvea2[4,1]
  PredictedNormal <- CCurvec+ (CCurved-CCurvec)/(1+(exp(-CCurvea*(Temperature-CCurveb))))
  Correction<-t((100-((Median-PredictedNormal)/(Median)*100))/100)
  i <- 1
  NormBothCorrect = matrix(NA, nrow=TotalRowsPer, ncol=TotalColsPer)
  repeat{
    j<-1
    repeat{
      NormBothCorrect[i,j]<-Data_Norm_Omit[i,j]*Correction[,j]
      j <- j+1
      if (j > TotalColsPer){
        break
      }
    }
    i<-i+1
    if (i > TotalRowsPer){
      break
    }
  }
  return(NormBothCorrect)
}
