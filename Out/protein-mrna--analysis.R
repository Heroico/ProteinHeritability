hause = read.table("hause-mrna-protein-correlation.txt")

png("hause-mrna-protein-r-squared.png")
plot(hause$V3 * hause$V3)
dev.off()

png("hause-mrna-protein-covariance.png")
plot(hause$V4)
dev.off()

hause_ext = read.table("hause-mrna-protein-ext-correlation.txt")

png("hause-mrna-protein-ext-r-squared.png")
plot(hause_ext$V3 * hause_ext$V3)
dev.off()

wu = read.table("wu-mrna-protein-correlation.txt")

png("wu-mrna-protein-r-squared.png")
plot(wu$V3 * wu$V3)
dev.off()
