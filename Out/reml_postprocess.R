library(ggplot2)

plot_heritability <- function(file_prefix){
	file <- paste(file_prefix, ".csv", sep="")
	data <- read.csv(file)
	ordered_data <- data[order(data$V.G._to_Vp),]

	sel <- data.frame(Index = seq.int(nrow(ordered_data)), Heritability = ordered_data$V.G._to_Vp)
	sel$y_min <- ordered_data$V.G._to_Vp - 2 * ordered_data$SE.V.G._to_Vp.
	sel$y_max <- ordered_data$V.G._to_Vp + 2 * ordered_data$SE.V.G._to_Vp.
	sel$nice <- ifelse(sel$y_min > 0,"yes", "no")

	p1<-ggplot(sel,aes(x=Index,y=Heritability, ymin = y_min, ymax=y_max) ) + 
				geom_pointrange(col='gray',alpha=0.7)+
				geom_point(aes(colour=nice))+
				coord_cartesian(ylim = c(-0, 1))+
				scale_colour_manual(values=c("yes" = "#991111", "no" = "#000000"))
	
	image <- paste(file_prefix, ".png", sep="")
	png(filename=image,width=1024,height=768)
	print(p1)
	dev.off()

	return()
}

plot_heritability("reml_results_hause")
plot_heritability("reml_results_wu")
