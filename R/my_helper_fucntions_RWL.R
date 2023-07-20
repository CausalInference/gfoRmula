
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





my.lag<-function(var.info){

  var.in.name <- var.info$var.name
  max.lag <- var.info$mx.lag
  my.output.list <-c()

  if (max.lag > 1){
    for(i in max.lag:2){
      my.output.list <- c(my.output.list, (paste("lag",i,"_",var.in.name,":=","lag",i-1,"_",var.in.name,sep="")))
    }
  }
  my.output.list<-c(my.output.list, (paste("lag1_",var.in.name,":=",var.in.name, sep="") ))

  return(my.output.list)
}

my.lag.function<-function(var.index){

  var.info <- lag_indicators[[var.index]]
  var.in.name <- var.info$var.name
  max.lag <- as.numeric(var.info$max.lag)


  if(max.lag > 0){
    if (max.lag > 1){
      for(i in max.lag:2){
        lhs <- paste0("lag",i,"_",var.in.name)
        rhs <- paste0("lag",i-1,"_",var.in.name)
        my.dt[,(lhs) := get(rhs)]
      }
    }
    lhs <- paste0("lag1_",var.in.name)
    my.dt[ , (lhs) := get(var.in.name)]
  }

}

my.lag.function2<-function(var.index){

  var.info <- lag_indicators[[var.index]]
  var.in.name <- var.info$var.name
  max.lag <- as.numeric(var.info$max.lag)


  if(max.lag > 0){
    if (max.lag > 1){
      for(i in max.lag:2){
        lhs <- paste0("lag",i,"_",var.in.name)
        rhs <- paste0("lag",i-1,"_",var.in.name)
        my.dt[,(lhs) := get(rhs)]
      }
    }
    lhs <- paste0("lag1_",var.in.name)
    my.dt[ , (lhs) := get(var.in.name)]
  }

}




my.lag.function3<- function(var.index,newdf=newdf){

  var.info <- lag_indicators[[var.index]]
  var.in.name <- var.info$var.name
  max.lag <- as.numeric(var.info$max.lag)
  if(max.lag > 0){
    max.lag2 <- max.lag - 1
    lag.names <- sapply(max.lag:1 , FUN=function(i){paste0("lag",i,"_",var.in.name)})
##    print(lag.names)

    newdf[,(lag.names) := lapply(max.lag2:0,FUN=function(x){
      if(x==0){eval(parse(text=var.in.name))}
      else {eval(parse(text=paste0("lag",x,"_",var.in.name)))}
    })]
  }
}


