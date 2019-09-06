getwd()
MyData <- read.csv(file="Desktop/Math/Cancer/CancerDataProject/alldat.csv", header=TRUE, sep=",")
#Curious about difference between RNA total count vs RNA observed count
#RNA_Ratio <- MyData[,4]/MyData[,9]
#RNA_Min <- min(RNA_Ratio)
gene_num <- 100

pois <- rpois(gene_num, 10)
#barplot(pois, main="pois")
hist(pois)

expo <- rexp(gene_num, 0.4)
#barplot(expo, main="expo")
hist(expo)

nbin <- rnbinom(gene_num, gene_num/2, 0.3)
#barplot(nbin, main="nbin")
hist(nbin)

hist(MyData[1:300,9])

plot(MyData[,2], MyData[,9])

nms <- c("Gene Number", "Total Count")

#PATIENT BP
df_bp <- data.frame(c(1:2554), MyData[1:2554, 9])
names(df_bp) <- nms
hist(df_bp[,2])
#plot(df_bp[,1], df_bp[,2])

#PATIENT OB
df_ob <- data.frame(c(2555:3537), MyData[2555:3537, 9])
names(df_ob) <- nms
hist(df_ob[,2])
#plot(df_ob[,1], df_ob[,2])

#PATIENT P07
df_p07 <- data.frame(c(3538:3771), MyData[3538:3771, 9])
names(df_p07) <- nms
hist(df_p07[,2])
#plot(df_p07[,1], df_p07[,2])

#PATIENT P72
df_p72 <- data.frame(c(3772:4070), MyData[3772:4070, 9])
names(df_p72) <- nms
hist(df_p72[,2])
#plot(df_p72[,1], df_p72[,2])

#PATIENT P35
df_p35 <- data.frame(c(4071:4346), MyData[4071:4346, 9])
names(df_p35) <- nms
hist(df_p35[,2])
#plot(df_p35[,1], df_p35[,2])

#PATIENT P40
df_p40 <- data.frame(c(4347:4554), MyData[4347:4554, 9])
names(df_p40) <- nms
hist(df_p40[,2])
#plot(df_p40[,1], df_p40[,2])
