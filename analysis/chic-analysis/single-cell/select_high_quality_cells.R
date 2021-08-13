library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(patchwork)


#Load count table of experiment
definecells<-read.csv("path_to_counttable/merged_count_table_50000.csv")

#correct wrongly importat cellnames and replace NA
definecells<-unite(definecells, names, c(1:3), remove=TRUE, sep = "-")
row.names(definecells)<-definecells$names
definecells <- definecells[,-1]
definecells[is.na(definecells)] <- 0

colnames(definecells)<-gsub(".", "-", colnames(definecells), fixed = TRUE)


#Load AT fraction table
tatable<-read.csv("path_to_TA_table/ScCHICLigation_merged_tagged.bamTA_obs_per_cell.csv")

#correct wrongly importat cellnames and replace NA
colnames(tatable)<-gsub(".", "-", colnames(tatable), fixed = TRUE)
rownames(tatable)<-tatable$X
tatable<-tatable[,-1]
tatable[is.na(tatable)] <- 0

#set n to the minimum read number per cell. It will get its counts per cell moved to a level where the TA fraction will become very low
#set m to the minimum %TA reads per cell. Choose both cutoff from QC plots
n=1000
m=70
tatable[2,]<-ifelse(tatable[2,] < n,1000000,tatable[2,])
testlist<-(as.numeric(tatable[1,])/as.numeric(tatable[2,])*100)>m
fraction<-tatable[,testlist]

#select cell based on TA table
table2<-definecells[,colnames(definecells) %in% colnames(fraction)]


#sort cells based on bin sums

i<-length(colnames(table2))
definecellsnorm<-table2

#normalize data
for(j in 1:i){
  definecellsnorm[,j]=table2[,j]/sum(table2[,j])*100000
}

orderlist<-order(rowSums(definecellsnorm,dims = 1), decreasing = TRUE)

definecellsordered<-table2[orderlist,]

#calcualte enrichment scores
length2<-length(rownames(definecellsordered))
m<-round(length2*0.2)
n<-round(length2*0.8)

results<-data.frame()
temp<-data.frame()
for(j in 1:i){
  temp[1,j]<-sum(definecellsordered[,j])
  temp[2,j]<-sum(definecellsordered[1:m,j])
  temp[3,j]<-sum(definecellsordered[n:length2,j])+1
  results[1,j]<-as.numeric((temp[2,j]/temp[1,j])*100)
  results[2,j]<-as.numeric((temp[3,j]/temp[1,j])*100)
  results[3,j]<-sum(definecellsordered[,j])
  results[4,j]<-results[1,j]/results[2,j]
}
results<-t(results)
colnames(results)<-c("fraction_of_reads_in_top_bins","fraction_of_reads_in_bottom_bins","total_reads","enrichment")

#plot results and select threshold for fraction of reads in top 20% bins
plot(results[,1],results[,2],xlab = "fraction of reads in top 20% bins", ylab = "fraction of reads in bottom 20% bins")
plot(results[,1],log10(results[,3]),xlab = "fraction of reads in top 20% bins", ylab = "log10 reads per cell")
hist(results[,1], breaks = 100, xlab = "fraction of reads in top 20% bins")

#set t to the lower threshold of the fraction of reads that should fall in the top20% of bins
t<-50
rownames(results)<-colnames(definecellsordered)
good20<-results[(results[,1]>t),]
output<-as.data.frame(results[,1])
rownames(output)<-rownames(results)
output[,1] <- c("bad")

output[rownames(output) %in% rownames(good20),1]<-"good"

#write table to be used for bamExtractSamples.py
write.table(output, file='filename', quote=FALSE, sep='\t', col.names = FALSE)
