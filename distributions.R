MyData <- read.csv(file="c:/Users/Preston/Desktop/Math/Cancer/Code/R/alldat.csv", header=TRUE, sep=",")

#Curious about difference between RNA total count vs RNA observed count
#RNA_Ratio <- MyData[,4]/MyData[,9]
#RNA_Min <- min(RNA_Ratio)
gene_num <- 100

pois <- rpois(gene_num, 10)
#barplot(pois, main="pois")
hist(pois)

expo <- rexp(gene_num, 0.9)
#barplot(expo, main="expo")
hist(expo)

nbin <- rnbinom(gene_num, gene_num/2, 0.3)
#barplot(nbin, main="nbin")
hist(nbin)

barplot(MyData[1:300,9])

data_frame <- data.frame(MyData[,4], MyData[,5])
