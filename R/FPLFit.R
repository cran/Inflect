#' This function determines the 4 parameter log fit constants for each protein in the melt shift data sets
#' @param Data_Norm_Omit List of accession numbers with NA omitted
#' @param NormBothCorrect List of normalized values from the FPLFit Correction function
#' @param Condition The Condition or the Control depending on which data set is being analyzed
#' @param Temperature The temperatures from the heat treatment
#' @param NumberTemperatures  The number of temperatures in the heat treatment
#' @importFrom stats optim
#' @return normalized data for each protein to the Inflect program
FPLFit<-function(Data_Norm_Omit,NormBothCorrect,Condition,Temperature,NumberTemperatures){

  TotalCellsPer<-as.numeric(NROW(Data_Norm_Omit)*NCOL(Data_Norm_Omit))
  TotalRowsPer<-as.numeric(NROW(Data_Norm_Omit))
  TotalColsPer<-as.numeric(NCOL(Data_Norm_Omit))

  n <- 1
  DataParameters = matrix(NA, nrow=TotalRowsPer, ncol=5)
  colnames(DataParameters) <- c(paste("a-",Condition),paste("Tm-",Condition),paste("c- ",Condition),paste("d-",Condition),paste("R2-",Condition))
  repeat{

    NormBothCorrectVal<-NormBothCorrect[n,1:TotalColsPer]


    df <- data.frame(x=Temperature,y=NormBothCorrect[n,])
    minmelt <- function(data, param) {with(data, sum(((param[3]/10+((param[4]-param[3]/10)/(1+exp(-param[1]*(x-param[2]*50)))))-y)^2))
    }

    CurveFit <- optim(par = c(-1,1,1,1), fn = minmelt, data = df, method=c('L-BFGS-B'),lower=-2,upper=2)

    CCurvea1<-CurveFit$par
    CCurvea<-CCurvea1[1]
    CCurveb<-CCurvea1[2]*50
    CCurvec<-CCurvea1[3]/10
    CCurved<-CCurvea1[4]

    DataParameters[n,1]<-CCurvea
    DataParameters[n,2]<-CCurveb
    DataParameters[n,3]<-CCurvec
    DataParameters[n,4]<-CCurved

    z<-1
    ssTot<-0
    ssResid<-0
    repeat{
      ssTot <- ((NormBothCorrectVal[z]-mean(NormBothCorrectVal))^2)+ssTot
      ssResid <- ((NormBothCorrectVal[z]-(CCurvec+((CCurved-CCurvec)/(1+exp(-CCurvea*(Temperature[z]-CCurveb))))))^2)+ssResid
      z<-z+1
      if (z > NumberTemperatures){
        break
      }
    }

    r2 <- 1 - (ssResid/ssTot)

    DataParameters[n,5]<-r2

    n=n+1
    if (n > TotalRowsPer){
      break
    }
  }
  return(DataParameters)
}
