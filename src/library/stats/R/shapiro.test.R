#  File src/library/stats/R/shapiro.test.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

shapiro.test <- function(x)
{
    DNAME <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if(is.na(n) || n < 3L || n > 5000L)
	stop("sample size must be between 3 and 5000")
    rng <- x[n] - x[1L]
    if(rng == 0) stop("all 'x' values are identical")
    if(rng < 1e-10) x <- x/rng # rescale to avoid ifault=6 with single version.
    res <- .Call(C_SWilk, x)
    RVAL <- list(statistic = c(W = res[1]), p.value = res[2],
		 method = "Shapiro-Wilk normality test", data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}

#exapnd the shapiro.test by adding the plotting feature to test normality
normarlity_test <- function(data,alpha=0.05,plotting = TRUE,verbose=TRUE) {
  
 if(any(is.na(data))) warning('There are missing values in your data. NAs will be omited in the follows steps')
  
 if (is.infinite(data)) stop("there is Inf value")
  
 if(verbose) message("Test in progress...")

 if(plotting == TRUE) {
    
    par(mfrow=c(1,2))
    
    #Histogram and Desity Plot
    hist(data,freq = FALSE,main = "Histogram and Density pLot",col="lightgray",na.rm=TRUE)
    lines(density(data,na.rm=TRUE),lty=1,col="red",lwd=2)
    
    #Normal density plot for comparison
    x<-c(round(min(data,na.rm = TRUE)):round(max(data,na.rm = TRUE)))
    lines(x,dnorm(x,mean(data,na.rm = TRUE),sd(data,na.rm = TRUE)),col="blue",lwd=2,lty=2)
    legend("topright", legend=c("Density", "Normality"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
    
    #QQ plot
    qqnorm(data,main="QQ-plot",col="blue")
    qqline(data,lty=1,col="red",lwd=2)
  }
  
    
    #Shapiro-Wilk's method
    result <- shapiro.test(data)
    p_value <- result$p.value
    
    if(p_value>alpha){
      print(paste("You data is normal,because your p_value =",p_value,">",alpha))
      
    }else{
      print(paste("You data is not normal,because your p_value =",p_value,"<=",alpha))
  
    
  }
}
