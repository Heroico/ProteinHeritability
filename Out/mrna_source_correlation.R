
mrna_source_data <- function(file_1, file_2) {
	data_1 = read.csv(file_1)
	data_2 = read.csv(file_2)

	genes <- names(data_1)

	correlation = numeric(length(genes))
	samples = numeric(length(genes))
	pvalues = numeric(length(genes))
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
		} else {
			correlation[index] = NA
		}
		samples[index] = count
		index = index+1
	}
	results <- data.frame(Gene = genes, Correlation = correlation, P.Value = pvalues, Samples = samples)
}

mrna_source_correlation <- function(data) {
	correlation = data$Correlation

	ordered_data <- correlation[order(correlation)]
	png(filename="predi_affy_by_gene.png",width=720,height=960)
	plot(ordered_data)
	dev.off()
	result = t.test(correlation, mu=0, conf.level=0.95)
}

data = mrna_source_data("../IntermediateC/intersection_mrna_affy.csv", "../IntermediateC/intersection_mrna_predi.csv")
write.table(data, "predixcan_affymetrix.table")
result = mrna_source_correlation(data)
print(result)
