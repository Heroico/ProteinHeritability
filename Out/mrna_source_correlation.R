library(ggplot2)

mrna_source_data <- function(file_1, file_2) {
	data_1 = read.csv(file_1)
	data_2 = read.csv(file_2)

	genes <- names(data_1)

	correlation = numeric(length(genes))
	samples = numeric(length(genes))
	pvalues = numeric(length(genes))
	y_min = numeric(length(genes))
	y_max = numeric(length(genes))
	index = 1
	for( i in genes ) {
		x <- data_1[[i]]
		y <- data_2[[i]]
		OK <- complete.cases(x, y)
		count <- length(OK[OK == TRUE])
		if (count > 6) {
			res = cor.test(x, y, alternative="two.sided", method="pearson",conf.level=0.95)
			pvalues[index] = res$p.value
			correlation[index] = res$estimate
			y_min[index] = res$conf.int[1]
			y_max[index] = res$conf.int[2]
		} else {
			correlation[index] = NA
		}
		samples[index] = count
		index = index+1
	}
	results <- data.frame(Gene = genes, 
			Correlation = correlation, P.Value = pvalues,
			 Samples = samples, Ymin = y_min, Ymax = y_max)
}

mrna_source_correlation <- function(data) {
	ordered_data <- data[order(data$Correlation),]
	sel <- data.frame(Index = seq.int(nrow(ordered_data)), 
					Correlation = ordered_data$Correlation,
					y_min = ordered_data$Ymin,
					y_max = ordered_data$Ymax)
	sel$nice <- ifelse(sel$y_min > 0,"yes", "no")

	p1<-ggplot(sel,aes(x=Index,y=Correlation, ymin = y_min, ymax=y_max) ) + 
				geom_pointrange(col='gray',alpha=0.1)+
				geom_point(aes(colour=nice))+
				coord_cartesian(ylim = c(-1, 1))+
				scale_colour_manual(values=c("yes" = "#991111", "no" = "#000000"))

	png(filename="predi_affy_by_gene.png",width=1024,height=768)
	print(p1)
	dev.off()
	result = t.test(data$Correlation, mu=0, conf.level=0.95)
}

qqunif = 
function(file = NULL,p,BH=T,CI=T,...)
{
  nn = length(p)
  xx =  -log10((1:nn)/(nn+1))
	if (! is.null(file)) {
		png(filename=file,width=1024,height=768)
	}
  plot( xx,  -sort(log10(p)),
     xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
       cex.lab=1.4,mgp=c(2,1,0),
       ... )
  abline(0,1,col='gray')
  if(BH)
    {
      abline(-log10(0.05),1, col='red',lty=1)
      abline(-log10(0.10),1, col='orange',lty=2)
      abline(-log10(0.25),1, col='yellow',lty=3)
      legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
             col=c('red','orange','yellow'),lty=1:3, cex=1)
      abline(h=-log10(0.05/nn)) ## bonferroni
    }
  if(CI)
  {
    ## create the confidence intervals
    c95 <- rep(0,nn)
    c05 <- rep(0,nn)
    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)
    ## this portion was posted by anonymous on
    ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
    
    for(i in 1:nn)
    {
      c95[i] <- qbeta(0.95,i,nn-i+1)
      c05[i] <- qbeta(0.05,i,nn-i+1)
    }
 
    lines(xx,-log10(c95),col='gray')
    lines(xx,-log10(c05),col='gray')

		if (! is.null(file)) {
			dev.off()
		}
  }
}


pa_data = mrna_source_data("../IntermediateC/intersection_mrna_affy.csv", "../IntermediateC/intersection_mrna_predi.csv")
write.table(pa_data, "predixcan_affymetrix.table")
result = mrna_source_correlation(pa_data)
print(result)
qqunif("qqunif_predixcan_affymetrix.png", pa_data$P.Value)


