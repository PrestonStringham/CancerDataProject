# Read in the data, and keep the names
BPO35T <- read.table("BPO35T.annot.top.Pathways.txt",sep="\t",header=TRUE,row.names=NULL)
OB9YER <- read.table("OB9YER.annot.top.Pathways.txt",sep="\t",header=TRUE,row.names=NULL)
P15007 <- read.table("P15007.annot.top.Pathways.txt",sep="\t",header=TRUE,row.names=NULL)
P15272 <- read.table("P15272.annot.top.Pathways.txt",sep="\t",header=TRUE,row.names=NULL)
P15335 <- read.table("P15335.annot.top.Pathways.txt",sep="\t",header=TRUE,row.names=NULL)
P15340 <- read.table("P15340.annot.top.Pathways.txt",sep="\t",header=TRUE,row.names=NULL)

dim(BPO35T) ## 2554   73
dim(OB9YER) ## 983  75
dim(P15007) ## 234  75
dim(P15272) ## 299  75
dim(P15335) ## 276  75
dim(P15340) ## 208  75

match(names(BPO35T),names(OB9YER))  #There are two missing
names(OB9YER)[is.na(match(names(OB9YER),names(BPO35T)))]
## "CHUANG_OXIDATIVE_STRESS_RESPONSE_DN" "ROETH_TERT_TARGETS_UP"              

diff(match(names(P15007),names(OB9YER)))
diff(match(names(P15272),names(OB9YER)))
diff(match(names(P15335),names(OB9YER)))
diff(match(names(P15340),names(OB9YER)))

# Let's clean up except for the columns we really want

shortnms <- c("BURTON_ADIPOGENESIS_3",
              "HALLMARK_MYC_TARGETS_V2",
              "HALLMARK_GLYCOLYSIS",
              "MOSERLE_IFNA_RESPONSE",
              "REACTOME_PYRIMIDINE_METABOLISM")

# nCount_RNA must be number of RNAs, which has information:
plot(BURTON_ADIPOGENESIS_3 ~ nCount_RNA,BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ nCount_RNA,BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ nCount_RNA,BPO35T))
## (Intercept) 4.12e-02   2.91e-03    14.2   <2e-16 ***
## nCount_RNA  1.33e-05   8.71e-07    15.2   <2e-16 ***

## "nCount_RNA"
## "total_counts" 
## "percent.mt" 
## "S.Score"  
## "G2M.Score" 
## seurat_clusters

plot(BURTON_ADIPOGENESIS_3 ~ total_counts,BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ total_counts,BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ total_counts,BPO35T))
## (Intercept)  5.27e-02   2.15e-03    24.6   <2e-16 ***
## total_counts 2.20e-06   1.33e-07    16.5   <2e-16 ***

plot(BURTON_ADIPOGENESIS_3 ~ log(total_counts),BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ log(total_counts),BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ log(total_counts),BPO35T))

plot(BURTON_ADIPOGENESIS_3 ~ log(nFeature_RNA),BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ log(nFeature_RNA),BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ log(nFeature_RNA),BPO35T))

plot(BURTON_ADIPOGENESIS_3 ~ percent.mt,BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ percent.mt,BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ percent.mt,BPO35T))
## (Intercept) 0.070998   0.002663   26.66  < 2e-16 ***
## percent.mt  0.000884   0.000214    4.13  3.7e-05 ***


plot(BURTON_ADIPOGENESIS_3 ~ S.Score,BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ S.Score,BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ S.Score,BPO35T))
## (Intercept) 0.081009   0.000959    84.4   <2e-16 ***
## S.Score     0.180582   0.003250    55.6   <2e-16 ***
# Highly non-linear and non-normal

plot(BURTON_ADIPOGENESIS_3 ~ G2M.Score,BPO35T)
abline(lm(BURTON_ADIPOGENESIS_3 ~ G2M.Score,BPO35T))
summary(lm(BURTON_ADIPOGENESIS_3 ~ G2M.Score,BPO35T))
## (Intercept)  0.08999    0.00105    85.9   <2e-16 ***
## G2M.Score    0.06837    0.00141    48.5   <2e-16 ***

plot(BURTON_ADIPOGENESIS_3 ~ as.factor(seurat_clusters),BPO35T)
summary(lm(BURTON_ADIPOGENESIS_3 ~ as.factor(seurat_clusters),BPO35T))
# Cluster 5 is a lot higher, for example

# Keep as correction factors
goodnms <-  c("nCount_RNA","nFeature_RNA","Sample","Patient","Time",
              "total_counts","percent.mt","S.Score","G2M.Score",
              "seurat_clusters","Encode_main_type","HMM_inferk4")

BPO35T.short <- BPO35T[,c(goodnms,shortnms)]
OB9YER.short <- OB9YER[,c(goodnms,shortnms)]
P15007.short <- P15007[,c(goodnms,shortnms)]
P15272.short <- P15272[,c(goodnms,shortnms)]
P15335.short <- P15335[,c(goodnms,shortnms)]
P15340.short <- P15340[,c(goodnms,shortnms)]

alldat <- rbind(
             cbind(patient="BP",clone="a",BPO35T.short),
             cbind(patient="OB",clone="a",OB9YER.short),
             cbind(patient="P07",clone="a",P15007.short),
             cbind(patient="P72",clone="a",P15272.short),
             cbind(patient="P35",clone="a",P15335.short),
             cbind(patient="P40",clone="a",P15340.short))
alldat$clone <- paste0(alldat$patient,alldat$HMM_inferk4)
names(alldat)[names(alldat)=="Patient"] <- "time"
patients <- unique(alldat$patient)

plot(BURTON_ADIPOGENESIS_3 ~ log(nFeature_RNA),alldat,pch=19,
     col=as.factor(patient),subset=nFeature_RNA > 2000)
for (ip in 1:length(patients)) {
  abline(lm(BURTON_ADIPOGENESIS_3 ~ log(nFeature_RNA),alldat,
         subset=patient==patients[ip] & nFeature_RNA > 2000),col=ip,lwd=3)
}
legend("topleft",patients,col=1:length(patients),lwd=3)
summary(lm(BURTON_ADIPOGENESIS_3 ~ log(nFeature_RNA)*as.factor(patient),alldat,
        subset=nFeature_RNA > 2000))

# Patients differ a lot in their counts
plot(log(nCount_RNA) ~ as.factor(patient),alldat)
plot(nFeature_RNA ~ as.factor(patient),alldat)
plot(log(total_counts) ~ as.factor(patient),alldat)

for (i in 1:length(shortnms)) {
  plot(alldat[,shortnms[i]] ~ as.factor(Time),alldat)
  title(main=shortnms[i])
  readline('hit return for next plot> ')
}

plot(HALLMARK_MYC_TARGETS_V2 ~ as.factor(clone),alldat)

myc.agg <- aggregate(HALLMARK_MYC_TARGETS_V2 ~ Time + clone,alldat,mean)

plot(myc.agg$HALLMARK_MYC_TARGETS_V2)

plot(nFeature_RNA ~ as.factor(Time),alldat)
summary(lm(nFeature_RNA ~ as.factor(Time),alldat))

plot(nFeature_RNA ~ as.factor(patient),alldat)
summary(lm(nFeature_RNA ~ as.factor(patient),alldat))

# Write it out
write.csv(alldat,file="alldat.csv")
