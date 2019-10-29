getwd()
alldat <- read.csv(file="Desktop/Math/Cancer/CancerDataProject/alldat.csv", header=TRUE, sep=",")

options(stringsAsFactors=FALSE)


#Make k cells
k <- 3

#m Genes
m <- 7

sim_start <- 1
sim_end <- k

#Create list for creating df_data
datalist = list()

#Total number of transcripts for each cell
numtrans <- c()

while(sim_start <= sim_end){
  #Select N uniformly on min and max of total counts
  N <- floor(runif(1, min=min(alldat[,4]), max=max(alldat[,4])))
  
  
  df_cells <- data.frame(transcript=paste0("t",1:m), q=runif(m))
  df_cells$q <- df_cells$q/sum(df_cells$q)
  
  tmp <- with(df_cells,sample(transcript,N,prob=q,replace=TRUE))
  tmp <- unlist(table(tmp))
  tmp <- data.frame(transcript=names(tmp),Ni=tmp)
  
  df_cells <- merge(df_cells,tmp[,c("transcript","Ni.Freq")],by="transcript",all.x=TRUE)
  df_cells$Ni.Freq[is.na(df_cells$Ni.Freq)] <- 0
  
  numtrans <- c(numtrans, sum(df_cells[,3]))
  
  p <- (1-df_cells[,2])^N
  var_F <- sum(p*(1-p))
  exp_F <- sum(1-(1-df_cells[,2])^N)
  
  dat <- df_cells
  datalist[[sim_start]] <- dat
  
  sim_start <- sim_start + 1
}

#Frame with all cells
df_data <-do.call(cbind, datalist)

#Auxiliary frame for computing GSEA scores
GSEAnms <- c("GENE_ID", "Transcripts", "Selected", "Rank")
df_GSEA <- data.frame(paste0("GENE",1:m), df_data[, 3], sample(c(0,1), replace=TRUE, size=m), rank(df_data[,3], ties.method="first"))
names(df_GSEA) <- GSEAnms

df_GSEA_rank_sorted <- df_GSEA[order(df_GSEA[,4]),]

alpha <- 0.1

i <- 1
e <- m
P_G_num <- c(matrix(0, 1, m))
P_NG_mult <- c(matrix(0, 1, m))
ng_count <- 0
while(i <= m){
  P_G_num[[i]] <- sum(df_GSEA_rank_sorted[1:i, 3]*df_GSEA_rank_sorted[1:i, 4])
  if(df_GSEA_rank_sorted[i, 3] == 0){
    ng_count <- ng_count + 1
    P_NG_mult[[i]] <- ng_count
  }
  P_NG_mult[[i]] <- ng_count
  i <- i + 1
}
df_G <- df_GSEA[df_GSEA[,3]==1,]
P_G <- (P_G_num)^alpha/(sum(df_G[,4]))^alpha
P_NG <- P_NG_mult*(1/(m-length(df_G[,1])))
ES <- sum(P_G - P_NG)

i <- 0
alpha_test <- seq(0, 1, by=0.1)
ES_list <- c()
while(i <= 10){
  ES_list <- c(ES_list, sum((P_G_num)^alpha_test[i]/(sum(df_G[,3]))^alpha_test[i] - P_NG))
  i <- i + 1
}
#plot(alpha_test, ES_list)
























#ADLER CODE
# This makes one cell
# df_cells <- data.frame(transcript=paste0("t",1:m), q=runif(m))
# df_cells$q <- df_cells$q/sum(df_cells$q)
# tmp <- with(df_cells,sample(transcript,N,prob=q,replace=TRUE))
# tmp <- unlist(table(tmp))
# tmp <- data.frame(transcript=names(tmp),Ni=tmp)
# print(dim(tmp))
# df_cells <- merge(df_cells,tmp[,c("transcript","Ni.Freq")],by="transcript",all.x=TRUE)
# df_cells$Ni.Freq[is.na(df_cells$Ni.Freq)] <- 0


# Scale up to make a whole bunch of cells where N can also be a random variable
# and for each cell extract summmary statistics such as N and F






# #With LP Scores
# pathways <- 3
# path_nms <- c(paste0("LP",1:pathways))
# while(sim_start < sim_end){
#   N <- floor(runif(1, min=min(alldat[,4]), max=max(alldat[,4])))
#   df_cells <- data.frame(transcript=paste0("t",1:m), q=runif(m))
#   df_cells$q <- df_cells$q/sum(df_cells$q)
#   
#   #tmp <- with(df_cells,sample(transcript,N,prob=q,replace=TRUE))
#   #tmp <- unlist(table(tmp))
#   #tmp <- data.frame(transcript=names(tmp),Ni=tmp)
#   #print(dim(tmp))
#   df_pathways <- data.frame(matrix(0, nrow = m, ncol = pathways))
#   names(df_pathways) <- path_nms
#   df_cells$path_nms[1:pathways] <- df_pathways[,1:pathways]
#   
#   
#   df_cells <- merge(df_cells,tmp[,c("transcript","Ni.Freq")],by="transcript",all.x=TRUE)
#   df_cells$Ni.Freq[is.na(df_cells$Ni.Freq)] <- 0
#   
#   p <- (1-df_cells[,2])^N
#   var_F <- sum(p*(1-p))
#   exp_F <- sum(1-(1-df_cells[,2])^N)
#   
#   
#   
#   sim_start <- sim_start + 1
# }