aic <- function( object, ... ){
  list0 <- list(object, ...)
  val0  <- c()
  for(i in 1:length(list0)){
    H    <- list0[[i]]
    if( ("sample.size" %in% names(H)) & ("n" %in% names(H)) )
        n <- H$sample.size
    if( ("sample.size" %in% names(H)) & !("n" %in% names(H)) )
        n <- H$sample.size
    if( !("sample.size" %in% names(H)) & ("n" %in% names(H)) )
        n <- H$n
    RSS  <- H$RSS
    p    <- length( H$par )
    temp <- 2*(p+1) + n*(log(2*pi)+1-log(n)+log(RSS))
    val0 <- c(val0, temp)  
  }
  return(val0)
}
