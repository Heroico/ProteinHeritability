mydata = read.table("mydata.txt")
v1 = mydata$V3
v2 = mydata$V4

png("pca.png")
plot(v1,v2)
dev.off()
