getwd()
alldat <- read.csv(file="Desktop/Math/Cancer/CancerDataProject/alldat.csv", header=TRUE, sep=",")

m <- 20000
df_cells <- data.frame(transcript=paste0("t",1:m), q=runif(m), Ni=0)


cell1 <- c(1:m)
cell2 <- c(1:m)
q <- c(matrix(1/m, 1, m))
detected <- c(matrix(0, 1, m))

df_cells <- data.frame(cell1, cell2, q, detected)

sample_start <- 1
N <- 1000
while(sample_start < N){
  
  select1 <- as.numeric(sample(1:m, 1))
  select2 <- as.numeric(sample(1:m, 1))
  if(!(select1==select2)){
    df_cells[select1, 4] <- df_cells[select1, 4] + 1
    df_cells[select2, 4] <- df_cells[select2, 4] + 1
  }
  
  sample_start <- sample_start + 1
  
}
exp_transcripts <- sum(mean(df_cells[1:m, 4]))
var_transcripts <- sum(var(df_cells[1:m, 4]))

hist(df_cells[,4])
hist(alldat[,5])
