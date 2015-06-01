library(ggplot2)

mrna_protein_correlation_results <- function(file_prefix) {
	mrna_file_name <- paste(file_prefix,"mrna_view.csv", sep="")
	mrna_data = read.csv(mrna_file_name)

	pheno_file_name <- paste(file_prefix,"pheno_view.csv", sep="")
	pheno_data = read.csv(pheno_file_name)

	genes <- names(mrna_data)
	genes <- genes[genes != "person_id"]
	
	data_len = length(genes)
	correlation = numeric(data_len)
	pvalue = numeric(data_len)
	ymin = numeric(data_len)
	ymax = numeric(data_len)
	index = 1	
	for (i in genes) {
		OK <- complete.cases(mrna_data[[i]], pheno_data[[i]])
		count <- length(OK[OK == TRUE])
		if (count > 6) {
			res = cor.test(mrna_data[[i]], pheno_data[[i]], alternative="two.sided", method="pearson",conf.level=0.95)
			correlation[index] = res$estimate
			pvalue[index] = res$p.value
			ymin[index] = res$conf.int[1]
			ymax[index] = res$conf.int[2] 
		} else {
			correlation[index] = NA
			pvalue[index] = NA
			ymin[index] = NA
			ymax[index] = NA
		}
		index = index +1
	}

	results <- data.frame(Gene = genes, Correlation = correlation, P.Value = pvalue, Ymin = ymin, Ymax = ymax)
}

plot_mrna_protein_correlation <- function(file_prefix) {
	results <- mrna_protein_correlation_results(file_prefix)

	ordered_data <- results[order(results$Correlation),]
	sel <- data.frame(Index = seq.int(nrow(ordered_data)), 
					Correlation = ordered_data$Correlation,
					y_min = ordered_data$Ymin,
					y_max = ordered_data$Ymax)
	
	sel$nice <- ifelse(sel$y_min > 0,"yes", "no")

	p1<-ggplot(sel,aes(x=Index,y=Correlation, ymin = y_min, ymax=y_max) ) + 
				geom_pointrange(col='gray')+
				geom_point(aes(colour=nice))+
				coord_cartesian(ylim = c(-1, 1))+
				scale_colour_manual(values=c("yes" = "#991111", "no" = "#000000"))
	
	image <- paste(file_prefix, "mrna_protein.png", sep="")
	png(filename=image,width=1024,height=768)
	print(p1)
	dev.off()

	return()
}

#TODO: remove
plot_mrna_protein_correlation("hause_")
plot_mrna_protein_correlation("wu_")

