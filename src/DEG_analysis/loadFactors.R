
#' @name readFactor
#'
#' @param configFile
#' @import NOISeq
#'
#' @return myFactors
#' @author Alberto Rodriguez-Izquierdo, 2021



readFactor <- function(configFile){
  
  factors <- read.csv(configFile$samplesFactorFilename, sep=';')

  return(factors)
}



#' @name loadFactors
#' @param factorFile
#' @param configFile
#' @import NOISeq
#'
#' @return myFactors
#' @author Alberto Rodriguez-Izquierdo, 2021


loadFactors <- function(factorFile,configFile){


  for (x in configFile$name_samples_1){

    if (exists('myfactor1')==FALSE){

      myfactor1 <- filter(factorFile, sample == x)

    }else{

      myfactor_x <- filter(factorFile, sample == x)

      myfactor1 <- rbind(myfactor1, myfactor_x)

    }

  }

  for (y in configFile$name_samples_2){

    if (exists('myfactor2')==FALSE){

      myfactor2 <- filter(factorFile, sample == y)

    }else{

      myfactor_y <- filter(factorFile, sample == y)

      myfactor2 <- rbind(myfactor2, myfactor_y)

    }
  }


  myfactor_first  <- rbind(myfactor1,myfactor2)

  myfactor <- data.frame(myfactor_first, batch= c(rep(1,nrow(myfactor1)), rep(2,nrow(myfactor2))))

  return(myfactor)
}

