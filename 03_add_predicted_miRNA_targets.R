setwd("~/Analysis/12_20190626_miRNA_Array/")

## read differential expression result table
res<-read.delim("analysis/result/conventional/res.IPscr_IPinh.txt")

########## add cluster
load("analysis/result/conventional/kmeans/10_clusters/k.rda")

res$cluster.kmeans.10<-k$cluster[match(as.character(res$SYMBOL), names(k$cluster))]


#################### get predicted miRNA targets #################


## get the gene names
gene.symbols<-as.character(res$SYMBOL)
gene.list<-strsplit(gene.symbols, fixed = T, split = "/")

names(gene.list)<-as.character(res$PROBEID)

# hash table
df.genes<-data.frame(SYMBOL=unlist(gene.list), PROBEID=rep(names(gene.list),
                sapply(gene.list,length)))


## read data

## miRDB
miRDB<-read.csv2("info/miRNA_predicted_targets/miRDB.csv")
mx<-match(as.character(df.genes$SYMBOL),as.character(miRDB$Gene.Symbol))

df.genes$miRDB<-as.character(miRDB$Gene.Symbol)[mx]

## miRMap
miRMap<-read.csv2("info/miRNA_predicted_targets/miRMap.csv")
mx<-match(as.character(df.genes$SYMBOL),as.character(miRMap$Gene))

df.genes$miRMap<-as.character(miRMap$Gene)[mx]

## miRWalk2
miRWalk2<-read.csv2("info/miRNA_predicted_targets/miRWalk2.0.csv")
mx<-match(as.character(df.genes$SYMBOL),as.character(miRWalk2$genesymbol))

df.genes$miRWalk2<-as.character(miRWalk2$genesymbol)[mx]


## targetScan
targetScan<-read.csv2("info/miRNA_predicted_targets/targetScan.csv")
mx<-match(as.character(df.genes$SYMBOL),as.character(targetScan$Target.gene))

df.genes$targetScan<-as.character(targetScan$Target.gene)[mx]



####### add to res.table
#miRDB
list.df.genes<-split(df.genes$miRDB, df.genes$PROBEID)
idx.miRDB.target<-sapply(list.df.genes,function(x) sum(!is.na(x))>0)

probeIDs.miRDB<-names(list.df.genes)[idx.miRDB.target]

res$miRDB<-(as.character(res$PROBEID) %in% probeIDs.miRDB)

#miRMap
list.df.genes<-split(df.genes$miRMap, df.genes$PROBEID)
idx.miRMap.target<-sapply(list.df.genes,function(x) sum(!is.na(x))>0)

probeIDs.miRMap<-names(list.df.genes)[idx.miRMap.target]

res$miRMap<-(as.character(res$PROBEID) %in% probeIDs.miRMap)

#miRWalk2
list.df.genes<-split(df.genes$miRWalk2, df.genes$PROBEID)
idx.miRWalk2.target<-sapply(list.df.genes,function(x) sum(!is.na(x))>0)

probeIDs.miRWalk2<-names(list.df.genes)[idx.miRWalk2.target]

res$miRWalk2<-(as.character(res$PROBEID) %in% probeIDs.miRWalk2)



#targetScan
list.df.genes<-split(df.genes$targetScan, df.genes$PROBEID)
idx.targetScan.target<-sapply(list.df.genes,function(x) sum(!is.na(x))>0)

probeIDs.targetScan<-names(list.df.genes)[idx.targetScan.target]

res$targetScan<-(as.character(res$PROBEID) %in% probeIDs.targetScan)


write.table(res, file="analysis/result/conventional/res.IPscr_IPinh_extended.txt",
            sep="\t", quote=F, row.names = F)


### get best cluster

cluster.split<-split(res, res$cluster.kmeans.10)


score.miRDB<-sapply(cluster.split, function(x) sum(x$miRDB)/nrow(x))
score.miRMap<-sapply(cluster.split, function(x) sum(x$miRMap)/nrow(x))
score.targetScan<-sapply(cluster.split, function(x) sum(x$targetScan)/nrow(x))
score.miRWalk2<-sapply(cluster.split, function(x) sum(x$miRWalk2)/nrow(x))

df.plot<-data.frame(data.base=rep(c("miRDB","miRMap","targetScan", "miRWalk2"),
                                  each=10),
        score=c(score.miRDB,score.miRMap, score.targetScan,score.miRWalk2),
        cluster=factor(rep(1:10,4)))

library(ggplot2)

pdf("analysis/result/conventional/miRNA_predicted_targets.pdf", height = 4,
    width=5)
    ggplot(df.plot, aes(cluster,score,fill=data.base))+
        geom_bar(stat = "identity", position = "dodge")+ 
        theme_bw()+ylab("fraction of members positive predicted")
dev.off()
