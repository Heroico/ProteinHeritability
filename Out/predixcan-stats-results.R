predixcan_results = read.table("predixcan-stats-results.txt")

png("predixcan-stats-results-r-squared.png")
plot(predixcan_results$V3 * predixcan_results$V3)
dev.off()
