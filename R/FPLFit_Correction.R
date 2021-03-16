#' This function determines the 4 parameter log fit constants for the two conditions in a replicate
#' @param Median The median fold change abundance from each temperature
#' @param Data_Norm_Omit List of accession numbers with NA omitted
#' @param Condition The Condition or the Control depending on which data set is being analyzed
#' @param Temperature The temperatures from the heat treatment
#' @importFrom stats nls
#' @return normalized data to the Inflect program
FPLFit_Correction<-function(Median,Data_Norm_Omit,Condition,Temperature){

  TotalCellsPer<-as.numeric(NROW(Data_Norm_Omit)*NCOL(Data_Norm_Omit))
  TotalRowsPer<-as.numeric(NROW(Data_Norm_Omit))
  TotalColsPer<-as.numeric(NCOL(Data_Norm_Omit))

  df <- data.frame(Temperature, Median)

  minmelt <- function(data, param) {with(data, sum(((param[3]+((param[4]-param[3])/(1+exp(-param[1]*(data[1]-param[2])))))-data[2])^2))
  }

  CurveFitCorrection <- optim(par = c(-1,mean(Temperature),min(Median),max(Median)), fn = minmelt, data = df, method=c('L-BFGS-B'))

  CCurvea1<-CurveFitCorrection$par
  CCurvea<-CCurvea1[1]
  CCurveb<-CCurvea1[2]
  CCurvec<-CCurvea1[3]
  CCurved<-CCurvea1[4]

  PredictedNormal <- CCurvec+ (CCurved-CCurvec)/(1+(exp(-CCurvea*(Temperature-CCurveb))))
  Correction<-((100-((Median-PredictedNormal)/(Median)*100))/100)

  i <- 1
  NormBothCorrect = matrix(NA, nrow=TotalRowsPer, ncol=TotalColsPer)
  repeat{
    j<-1
    repeat{
      NormBothCorrect[i,j]<-Data_Norm_Omit[i,j]*Correction[j]
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
