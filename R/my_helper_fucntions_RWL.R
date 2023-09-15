
update_lag_indicator_rwl <- function(models, covs,  lag_indicator){
  n<- length(covs)


  test1 <- (stringr::str_extract_all(deparse(models),'lag\\d[_]\\w+'))

  yy <- (unique(unlist(test1)))
  ##print(yy)

  ##print(covs)
  mytest <- sapply((covs),function(var){

    xx <- all(is.na( stringr::str_match(string=yy,
                                        pattern = paste("_", var,sep=""))))

    if (xx == FALSE){
      z1<-(unlist(stringr::str_extract_all(
        string=yy,
        pattern = paste("lag\\d+[_]",var,
                        sep=""))))

      z1.lags <- stringr::str_extract_all(string=z1 ,
                                          pattern ='lag\\d+')
      z1.lags <- as.numeric(stringr::str_extract_all(string = z1.lags ,
                                                     pattern = '\\d+'))
      z1.lags <- max(z1.lags)
    } else {
      z1.lags <- 0
    }
    return(list(var.name = var, max.lag=z1.lags))
    #return(z1.lags)

  } , simplify = FALSE , USE.NAMES = TRUE)


  return(mytest)
}

update_lagcumavg_indicator_rwl <- function(models, covs,  lag_indicator){
  n<- length(covs)

  test1 <- (stringr::str_extract_all(deparse(models),'lag_cumavg\\d[_]\\w+'))

  yy <- (unique(unlist(test1)))
  mytest <- sapply((covs),function(var){

    xx <- all(is.na( stringr::str_match(string=yy,
                                        pattern = paste("_", var,sep=""))))

    if (xx == FALSE){
      z1<-(unlist(stringr::str_extract_all(
        string=yy,
        pattern = paste("lag_cumavg\\d+[_]",var,
                        sep=""))))

      z1.lags <- stringr::str_extract_all(string=z1 ,
                                          pattern ='lag_cumavg\\d+')
      z1.lags <- as.numeric(stringr::str_extract_all(string = z1.lags ,
                                                     pattern = '\\d+'))
    } else {
      z1.lags <- 0
    }


    needed.list <-NULL

    needed.list <- sapply(z1.lags , FUN=function(i){if(i > 1 ){ paste0("lag",i-1,"_",var)}})

    return(list(var.name = var, max.lag=z1.lags, needed.list = needed.list))


  } , simplify = FALSE , USE.NAMES = TRUE)

  return(mytest)
}


update_cumavg_indicator_rwl <- function(models, covs,  lag_indicator){
  n<- length(covs)

  test1 <- (stringr::str_extract_all(deparse(models),'cumavg_\\w+'))

  yy <- (unique(unlist(test1)))

  mytest <- sapply((covs),function(var){

    xx <- all(is.na( stringr::str_match(string=yy,
                                        pattern = paste("cumavg_", var,sep=""))))

    if (xx == FALSE){
      z1.lags <- 1
    } else {
      z1.lags <- -1
    }
    return(list(var.name = var, max.lag=z1.lags))

  } , simplify = FALSE , USE.NAMES = TRUE)


  return(mytest)
}


update.lagged.for.lag_cumavg.v2<-function(lag_indicators,needed.list,covnames){

  mylag.test<-lapply(seq_along(covnames ),FUN=function(index){
    # there should be a 1-1 mapping of the ordering of covnames and lag_indicators

    tmp.lag.list <- lag_indicators[[index]]
    lagged.max.lag <- tmp.lag.list$max.lag

    for(i in seq_along(needed.list)){

      lag.check <- needed.list[i]
      lag.check.name <- unlist(strsplit(lag.check,"_"))[2]
      lag.part <- unlist(strsplit(lag.check,"_"))[1]
      needed.lag <- as.numeric((stringr::str_extract_all(string=lag.part,pattern='\\d')))

      if (tmp.lag.list$var.name == lag.check.name && tmp.lag.list$max.lag < needed.lag) {

        tmp.lag.list$max.lag <- needed.lag
      }
    }
    lag_indicators[[index]] <- tmp.lag.list
  })

  return(mylag.test)
}

find.max.lag.needed <- function(lag_indicators){
  max.lag <- 0
  for (i in seq_along(lag_indicators)){
    if (lag_indicators[[i]]$max.lag > max.lag) {max.lag <- lag_indicators[[i]]$max.lag}
  }

  return(max.lag)

}



