
# Changed space to underline in this file
genesets <- read.table("AOCS093_10X_pathways_genesets.txt",
                            header=TRUE,row.names=NULL)
genesets <- genesets[,1:3]
max(genesets$Total_Genes)
# 1972 seems a bit large, but such is life
# This makes a good check with our data

# This came from stackoverflow in a big way.  
inlist <- strsplit(readLines("c2.all.v7.0.symbols.gmt"), "[[:space:]]+")
c2.pathways <- lapply(inlist, tail, n = -2)
names(c2.pathways) <- lapply(inlist, head, n = 1)

inlist <- strsplit(readLines("c6.all.v7.0.symbols.gmt"), "[[:space:]]+")
c6.pathways <- lapply(inlist, tail, n = -2)
names(c6.pathways) <- lapply(inlist, head, n = 1)

inlist <- strsplit(readLines("h.all.v7.0.symbols.gmt"), "[[:space:]]+")
h.pathways <- lapply(inlist, tail, n = -2)
names(h.pathways) <- lapply(inlist, head, n = 1)

sum(is.na(match(names(c2.pathways),genesets$Gene_Set))) #1191 missing
sum(is.na(match(names(c6.pathways),genesets$Gene_Set))) #189 missing
sum(is.na(match(names(h.pathways),genesets$Gene_Set))) #1 missing

# Any repeats?
intersect(names(c2.pathways),names(c6.pathways))
intersect(names(c2.pathways),names(h.pathways))
intersect(names(h.pathways),names(c6.pathways))
# No repeats, glom together into one file

all.pathways <- c(c2.pathways,c6.pathways,h.pathways)
sum(is.na(match(names(all.pathways),genesets$Gene_Set))) #1381 missing

length(all.pathways)
# 5740
dim(genesets)
# 4774

# Not quite sure about this, but we'll do our best with the ones that
# match.

tmp <- lapply(all.pathways,length)
tmp <- unlist(tmp)

tmp <- data.frame(Gene_Set=names(tmp),gcount=tmp)
tmp <- merge(genesets,tmp,by="Gene_Set",all.x=TRUE,all.y=TRUE)
plot(gcount ~ Total_Genes,tmp,pch=19)
abline(0,1)
# These don't quite match for some reason.

# There are 1381 in all.pathways that aren't in our Gene_Sets, and 
# 415 in Gene_Sets that are not in all.pathways.

goodpathways <- tmp$Gene_Set[is.na(tmp$Total_Genes + tmp$gcount)==FALSE]
# There are 4359 in both

# Need to figure out how to pull these out of the list, and then make
# the list into a giant matrix

good.pathways <- all.pathways[goodpathways]
# OK, that was easy

tmp <- unlist(all.pathways)
tmp <- unique(tmp)
length(tmp)
allgenes <- tmp
# 20817.  That's a lot

tmp <- unlist(good.pathways[2])
match(tmp,allgenes)
tfun <- function(x) match(unlist(x),allgenes)
tmp <- lapply(good.pathways,tfun)
# This pulls out all the numbers

hugemat <- matrix(0,length(allgenes),length(good.pathways))
# All we need to do is fill in 1's in the rows denoted by tmp[i] in
# column i

# We could do this as a loop...
hugemat[as.numeric(unlist(tmp[1])),1] <- 1

for (i in 1:length(tmp)) {
  hugemat[as.numeric(unlist(tmp[i])),i] <- 1
}
# That was fast
# Make it into a data frame
hugemat <- as.data.frame(hugemat)
dimnames(hugemat)[[1]] <- allgenes
dimnames(hugemat)[[2]] <- names(good.pathways)

# We may have to subset hugemat for the genes we have in the datasets,
# but this should have all the information

tbl1 <- table(apply(hugemat,1,sum))
# 289 genes are in no pathway, 1 gene is in 260 pathways
plot(log(1+tbl1))
# Number of pathways per gene very close to exponential except for 0 and 1.
plot(log(1+tbl1),xlim=c(0,80))

tbl2 <- table(apply(hugemat,2,sum))
plot(tbl2)
plot(log(1+tbl2))
plot(log(1+tbl2),xlim=c(0,200))
# Strange shape. 
# Number of genes per pathway definitely not exponential

# Might as well use AOCS093 since the pathwyas are already here.
AOCS093.counts <- read.table("AOCS093.10X.counts.txt",
             header=TRUE,row.names=NULL)
# It read fine. No surprises there.
# Now time to compute scores.
dim(AOCS093.counts)
#  16687  4255

dim(hugemat)
# 20817  4359
# They differ but not by too much. Likely will have NA genes.

gene_match_patient <- match(AOCS093.counts[,1], row.names(hugemat))
length(gene_match_patient)
# 16687

gene_match_patient <- gene_match_patient[!is.na(gene_match_patient)]
length(gene_match_patient)
# 15700, lost about a thousand genes
table(table(gene_match_patient))
#    1 
# 15700    -   Nice, no repetitions. 

subset <- hugemat[gene_match_patient, ]
dim(subset)
# 15700  4359
# Missing columns

# These are in counts per million:
sum(contest_data4[,2])

# Deleted two rows with /// that looked weird
contest_data2 <- read.table("cleandataset2.count.txt",header=TRUE,row.names=NULL)

data2_genes <- contest_data2[,1]

data2_indices <- match(data2_genes, allgenes)
length(data2_indices)
# 22240

data2_indices <- data2_indices[!is.na(data2_indices)]
length(data2_indices)
# 18794

table(table(data2_indices))
#     1     2 
# 18790     2

data2_indices <- unique(data2_indices)
length(data2_indices)
# 18792

data2_goodgenes <- allgenes[data2_indices]

data2_matrix <- matrix(0, length(data2_goodgenes), length(all.pathways))

data2_scores <- matrix(0, length(contest_data2[1,])-1, length(all.pathways))

# I think this will work
for(i in c(1:length(all.pathways))){
	tmp_indices <- match(all.pathways[[i]], data2_goodgenes)
	tmp_indices <- tmp_indices[!is.na(tmp_indices)]
	
	data2_matrix[tmp_indices, i] <- 1
	score <- apply(contest_data2[tmp_indices,2:data2_length],2,sum)/apply(contest_data2[,2:data2_length],2,sum)
	data2_scores[,i] <- score
}
# this took 30 minutes to compute. Need to fond another way. 

data2_df <- as.data.frame(data2_scores)
colnames(data2_df) <- names(all.pathways)



data2_length <- ncol(contest_data2)
# data2_vector <- apply(contest_data2[data2_indices,2:data2_length],2,sum)/apply(contest_data2[,2:data2_length],2,sum)
data2_score_mat <- matrix(data2_vector, length(data2_vector), length(all.pathways))
data2_multiplied <- data2_matrix * data2_score_mat
data2_df <- as.data.frame(data2_multiplied)
rownames(data2_df) <- allgenes[data2_indices, 1]

contest_data4 <- read.table("dataset4.tpm.txt",
                           header=TRUE, row.names=NULL)
# Crossing fingers
# It didn't give me an error!

# No idea what these numbers mean. Checking email again. Checking didn't
# help. Interesting that numbers aren't integers. Let's read the second
# dataset in anyway.

data4_indices <- match(contest_data4[,1], allgenes)
#Anything we can salvage? Looks like some.

data4_indices <- data4_indices[!is.na(data4_indices)]

table(table(data4_indices))
#     1     2     3     4     5     6     7    11    13    30 
# 19137   117     3     5     3     3     1     1     1     1    -   yikes

#unique it, I guess?
data4_indices <- unique(data4_indices)

table(table(data4_indices))
#     1 
# 19272 

# I don't know if it is worth continuing. I just lost well over half of the genes from the original dataset. 
# Let's talk.

# Nov 22: We just need to figure out which genes are in our original
# genesets, and then add up the fractions.  The Counts Per Million are
# already in the naively corrected format, so we divide each column by
# 10^6, and multiply by hugemat.  I'm a bit worried that  contest_data4
# has 57823 rows, and looking at some of the names (such as the ones
# starting with RP11, there might be a bunch of repeats).

table(substr(contest_data4[,1],1,4)=="RP11")
##FALSE  TRUE 
##46099 11724 

table(substr(contest_data2[,1],1,4)=="RP11")
##FALSE 
##22240 


# There are weird
table(table(substr(contest_data4[,1],1,4)))
table(table(substr(contest_data2[,1],1,4)))

table(table(substr(contest_data4[,1],1,5)))
table(table(substr(contest_data2[,1],1,5)))

table(table(substr(dimnames(hugemat)[[1]],1,5)))

# I tried to look up what RP11  is. Looks like RP stands for Roswell
# Park. As in, Roswell Park Cancer Center.  I am still a little confused
# about what they are. Google search for individual RP11s resulted in
# databases Saying they are genes, but I also heard they could be
# lnc-RNA from samples taken at Roswell Park, which is how it got its
# name. For now, in the above sets, I just have been looking for exact
# matching names of genes, but 57823 is concerning. 

# This looks pretty good, but let's try it with a matrix multiply.

matchgenes2 <- intersect(contest_data2$Gene.Symbol,rownames(hugemat))
tmpdat <- subset(contest_data2,Gene.Symbol %in% matchgenes2)
# Grrr.  Duplicated rows in contest_data2
tmp <- subset(contest_data2,duplicated(Gene.Symbol))
badgenes <- tmp$Gene.Symbol
tmp <- subset(contest_data2,Gene.Symbol %in% badgenes)
apply(tmp[,2:ncol(tmp)],1,sum)
# Well, these don't match at all, but are kind of small.  We are going
# to just keep the first one.

cleandat2 <- subset(contest_data2,duplicated(Gene.Symbol)==FALSE)
tmpdat <- subset(cleandat2,Gene.Symbol %in% matchgenes2)
tmpmat <- as.matrix(hugemat[matchgenes2,])
tmpdat <- t(as.matrix(tmpdat[,2:ncol(tmpdat)]))
matdat2 <- tmpdat%*%tmpmat
# Forgot to normalize by dividing by the read count
dmat <- diag(1/apply(tmpdat,1,sum))
matdat2 <- dmat%*%matdat2
matdat2 <- as.data.frame(matdat2)
# Names vanished.  Need to restore
dimnames(matdat2)[[1]] <- names(cleandat2[,2:ncol(cleandat2)])

# Let's compare
dim(data2_df) 

pathnm <- "ABBUD_LIF_SIGNALING_1_DN"
plot(matdat2[,pathnm],data2_df[,pathnm],pch=19)
pathnms <- names(matdat2)
pathnm <- pathnms[234]
plot(matdat2[,pathnm],data2_df[,pathnm],pch=19)

# Darn.  They don't match.  Need to do one by hand
v1 <- tmpmat[,pathnm]
v2 <- tmpdat[1,]
sum(v1*v2)/sum(v2)
# Woohoo! This matches matdat2[1,pathnm]
# Boohoo. This does not match data2_df[1,pathnm]

# On to contest_data4

matchgenes4 <- intersect(contest_data4$gene_name,rownames(hugemat))
tmpdat <- subset(contest_data4,gene_name %in% matchgenes4)
# Grrr.  Duplicated rows in contest_data4
tmp <- subset(contest_data4,duplicated(gene_name))
badgenes <- tmp$gene_name
tmp <- tmp[order(tmp$gene_name),]
tmp <- subset(contest_data4,gene_name %in% badgenes)
tmp2 <- data.frame(gene_name=tmp$gene_name)
tmp2$sum <- unlist(as.vector(unlist(apply(tmp[,2:ncol(tmp)],1,sum))))
tmp2 <- tmp2[order(tmp2$gene_name),]

# Well, these don't match at all, and are all over the place, often
# including a single large one.  We are going to keep that.
tmp <- contest_data4
tmp$sum <- unlist(as.vector(unlist(apply(tmp[,2:ncol(tmp)],1,sum))))
tmp <- tmp[order(tmp$gene_name,-tmp$sum),]
cleandat4 <- subset(tmp,duplicated(gene_name)==FALSE)
cleandat4$sum <- NULL

tmpdat <- subset(cleandat4,gene_name %in% matchgenes4)
tmpmat <- as.matrix(hugemat[matchgenes4,])
tmpdat <- t(as.matrix(tmpdat[,2:ncol(tmpdat)]))
matdat4 <- tmpdat%*%tmpmat
# Forgot to normalize by dividing by the read count
dmat <- diag(1/apply(tmpdat,1,sum))
# These should all have been 1.0e-6, but are off because we eliminated
# the duplicated genes.
matdat4 <- dmat%*%matdat4
# Names vanished.  Need to restore
matdat4 <- as.data.frame(matdat4)
dimnames(matdat4)[[1]] <- names(cleandat4[,2:ncol(cleandat4)])

# Before writing these out to submit for our prize, let's graph
# something

par(mfrow=c(3,3))
prand <- sample(pathnms,9)
for (pathnm in prand) {
  xrange <- range(c(matdat2[,pathnm],matdat4[,pathnm]))
  plot(ecdf(matdat2[,pathnm]),xlim=xrange,lwd=2,main="")
  lines(ecdf(matdat4[,pathnm]),lwd=2,col="red")
  title(main=pathnm)
}
par(mfrow=c(1,1))
# These look rather interesting

# Let's write 'em out!
write.csv(matdat2,"stringham_adler_dataset2.csv")
write.csv(matdat4,"stringham_adler_dataset4.csv")

# Test with one pathway.
hugemat_small <- as.matrix(hugemat[,1])
rownames(hugemat_small) <- rownames(hugemat)
colnames(hugemat_small) <- colnames(hugemat[1,][1])

cleandat2_small <- subset(contest_data2,duplicated(Gene.Symbol)==FALSE)
tmpdat_small <- subset(cleandat2,Gene.Symbol %in% matchgenes2)
tmpmat_small <- as.matrix(hugemat_small[matchgenes2,])
tmpdat_small <- t(as.matrix(tmpdat_small[,2:ncol(tmpdat_small)]))
matdat2_small <- tmpdat_small%*%tmpmat_small
dmat_small <- diag(1/apply(tmpdat_small,1,sum))
matdat2_small <- dmat_small%*%matdat2_small
matdat2_small <- as.data.frame(matdat2_small)
dim(matdat2_small)
# [1] 626   1
# It worked!


# Let's try to make hugemat stored as a sparse matrix.
require(Matrix)
compressed_hugemat <- Matrix(as.matrix(hugemat), sparse=TRUE)

object.size(hugemat)
# 727947560 bytes

object.size(compressed_hugemat)
# 6991520 bytes

# Looks smaller
sys.time()
cleandat2 <- subset(contest_data2,duplicated(Gene.Symbol)==FALSE)
tmpdat <- subset(cleandat2,Gene.Symbol %in% matchgenes2)
tmpmat <- as.matrix(compressed_hugemat[matchgenes2,])
tmpdat <- t(as.matrix(tmpdat[,2:ncol(tmpdat)]))
matdat2 <- tmpdat%*%tmpmat
# Forgot to normalize by dividing by the read count
dmat <- diag(1/apply(tmpdat,1,sum))
matdat2 <- dmat%*%matdat2
matdat2 <- as.data.frame(matdat2)
