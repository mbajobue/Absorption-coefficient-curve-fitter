#!/usr/bin/env Rscript

# First, the data is imported
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
	data = read.csv(args[1]) 
} else {
	data = read.csv("data.csv")
}

z_max <- ncol(data)
z <- 1 # Counter that indicates which set of data is analysed

# Plot and output image file parameters
png( "plot.png", width = 5.55, height = 5.55, units = "in", res = 300, pointsize = 6)
par( mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)

plot(data[[1]],data[[2]], pch=4, type="n", ylab = expression(paste(alpha," (cm"^"-1",")")),xlab = expression(paste("Energy (eV)")))



# First estimation of eg
est_eg <- function(z){

	b <- data[[z+1]][!is.na(data[[z+1]])]
	b <- b[b>max(b)/2]
	a <- data[[z]][match(b,data[[z+1]])]


	fit <- lm(a ~ b)


	eg <- coef(fit)[[1]]
	return(eg)
}


### Main function ###


repeat{
	en <- 0.05 	# A small value of en is assumed
	eg <- est_eg(z)	
	ec <- (eg+2*en)	
	
	# alpha1 and alpha 2 starting values.
	alpha2 <- 100000
	alpha1 <- 10
	
	
	repeat{
	
		# NA values from our data are eliminated
		x <- data[[z]][!is.na(data[[z]])]
		y <- data[[z+1]][!is.na(data[[z+1]])]
		
		# For e>ec:
		x1 <- x[x>ec]
		y1 <- head(y, length(x1))
		F1 <- function(x1,a2,egg){a2*(x1-egg)^2}
		fit1 <- nls( y1 ~ F1(x1,a2,egg),start=list(a2=alpha2,egg=eg)) 
		
		# New estimation of eg and the relative difference with the previous estimation
		rel=abs(coef(fit1)[[2]]-eg)/coef(fit1)[[2]]
		eg <- coef(fit1)[[2]] 
		
		
		# For e<ec:
		x2 <- x[x<ec]  
		y2 <- tail(y, length(x2))
		F2 <- function(x2,a1){a1*(exp((x2-eg)/en))}
		fit2 <- nls( y2 ~ F2(x2,a1),start=list(a1=alpha1))
		
		
		# New estimated values
		alpha1 <- coef(fit2)[[1]]
		alpha2 <- coef(fit1)[[1]]
		en <- sqrt(alpha1/(4*alpha2)*exp(2)) 
		ec <- eg+2*en		
		
		
		
	if(rel < 10^(-5)){
		break
	  }
	}
	
	points(data[[z]],data[[z+1]], pch=4)                    
	lines(smooth.spline(c(x1,head(x2,n=1)),c(predict(fit1),
	head(predict(fit2),n=1))),lwd=1,col='red')              
	lines(smooth.spline(x2,predict(fit2)),lwd=1,col='red')  
	
	
	
	# Errors
	err_eg <- coef(summary(fit1))[,2][[2]]
	err_alpha2 <- coef(summary(fit1))[,2][[1]]
	err_alpha1 <- coef(summary(fit2))[,2][[1]]
	err_en <- (exp(1)/2)*sqrt((err_alpha1/alpha2)^2+((alpha1*err_alpha2)/alpha2^2)^2)
	err_ec <- sqrt(err_eg^2+4*(err_en)^2)
	
	# The eg and its error is printed out
	cat("Experiment",paste((z+1)/2),":\n","eg = ",eg,"±",signif(err_eg,1),"\n")


	# A csv file is generated with the estimated values and their errors
	est_val <- c(eg,en,ec,alpha1,alpha2)
	errors <- c(err_eg,err_en,err_ec,err_alpha1,err_alpha2) 
	df <- data.frame(est_val,errors)
	write.csv(df, file=print(paste("Experiment_",(z+1)/2,".csv",sep="")),row.names=c("eg","en","ec","alpha1","alpha2"))
	
	
	z <- z+2 # The counter z is set to analyse the next experiment

if(z>z_max){
	break
  }
}

dev.off()
