
mrna_source_correlation <- function(file_1, file_2) {
	data_1 = read.csv(file_1)
	data_2 = read.csv(file_2)

	genes <- names(data_1)

	correlation = numeric(length(genes))
	index = 1
	for( i in genes ) {
		x <- data_1[[i]]
		y <- data_2[[i]]
		OK <- complete.cases(x, y)
		count <- length(OK[OK == TRUE])
		if (count > 6) {
			res = cor.test(x, y, alternative="two.sided", method="pearson",conf.level=0.95)
			correlation[index] = res$estimate
		} else {
			correlation[index] = NA
		}
		index = index+1
	}

	ordered_data <- correlation[order(correlation)]
	png(filename="predi_affy_by_gene.png",width=720,height=960)
	plot(ordered_data)
	dev.off()
	result = t.test(correlation, mu=0, conf.level=0.95)
}

result = mrna_source_correlation("../IntermediateC/intersection_mrna_affy.csv", "../IntermediateC/intersection_mrna_predi.csv")
print(result)
