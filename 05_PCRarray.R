setwd("/Users/admin/Library/Mobile Documents/com~apple~CloudDocs/Analysis_sync/S017-miR132/")
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

array<-read.csv2("PCR Array-FC und Pvalue.csv")


######## SFCplus EXT ########
df<-array[,c(1,grep("SFCplus_Acq", colnames(array)))]
colnames(df)<-c("Gene","logFC","p.value")
df$logFC<-as.numeric(df$logFC)
df$p.value<-as.numeric(df$p.value)

df$pVal<-"ns"
df$pVal[df$p.value<0.07]<-"< 0.07"
df$pVal[df$p.value<0.05]<-"< 0.05"



g<-ggplot(df,aes(x=logFC,y=-log10(p.value), fill=pVal) )+
    geom_point(pch=21,size=2)+theme_classic()+
    scale_fill_manual(values = c("darkorange", "grey","white" ))+
    ggtitle("SFC+/ACQ")

#get top 20 p-adjusted value
df.top<-subset(df, pVal!="ns")


pdf("SFCplus_Acq_Volcano.pdf", width=4, height = 3.5)
    print(g+geom_text_repel(data = df.top,
                        aes(label = Gene),
                        size = 3,
                        box.padding = unit(0.3, "lines"),
                        point.padding = unit(0.2, "lines")))
dev.off()


######## SFCminuns EXT ########


df<-array[,c(1,grep("SFCminus_Ext", colnames(array)))]
colnames(df)<-c("Gene","logFC","p.value")
df$logFC<-as.numeric(df$logFC)
df$p.value<-as.numeric(df$p.value)

df$pVal<-"ns"
df$pVal[df$p.value<0.07]<-"< 0.07"
df$pVal[df$p.value<0.05]<-"< 0.05"


g<-ggplot(df,aes(x=logFC,y=-log10(p.value), fill=pVal) )+
    geom_point(pch=21,size=2)+theme_classic()+
    scale_fill_manual(values = c("darkorange", "grey","white" ))+
    ggtitle("SFC-/EXT")

#get top 20 p-adjusted value
df.top<-subset(df, pVal!="ns")


pdf("SFCminus_Ext_Volcano.pdf", width=4, height = 3.5)
    print(g+geom_text_repel(data = df.top,
                        aes(label = Gene),
                        size = 3,
                        box.padding = unit(0.3, "lines"),
                        point.padding = unit(0.2, "lines")))
dev.off()


######## SFCplus EXT ########


df<-array[,c(1,grep("SFCplus_Ext", colnames(array)))]
colnames(df)<-c("Gene","logFC","p.value")
df$logFC<-as.numeric(df$logFC)
df$p.value<-as.numeric(df$p.value)

df$pVal<-"ns"
df$pVal[df$p.value<0.07]<-"< 0.07"
df$pVal[df$p.value<0.05]<-"< 0.05"


g<-ggplot(df,aes(x=logFC,y=-log10(p.value), fill=pVal) )+
    geom_point(pch=21,size=2)+theme_classic()+
    scale_fill_manual(values = c("darkorange", "grey","white" ))+
    ggtitle("SFC+/EXT")

#get top 20 p-adjusted value
df.top<-subset(df, pVal!="ns")


pdf("SFCplus_Ext_Volcano.pdf", width=4, height = 3.5)
    print(g+geom_text_repel(data = df.top,
                        aes(label = Gene),
                        size = 3,
                        box.padding = unit(0.3, "lines"),
                        point.padding = unit(0.2, "lines")))
dev.off()


## heatmap
library(pheatmap)


mx<-as.matrix(array[,grep("Fold", colnames(array))])
mx<-apply(mx,2,as.numeric)
genes<-array$Gene
genes[duplicated(genes)]<-paste0(genes[duplicated(genes)], ".2")

rownames(mx)<-genes
colnames(mx)<-c("SCF+/Acq","SCF-/Ext","SCF+/Ext" )

color.heat3 = colorRampPalette(rev(brewer.pal(n = 7,
                                 name ="RdBu"))[c(rep(1,3),1:7,rep(7,3))])(100)

pdf("PCRarray_heatmap.pdf", height = 12, width = 4)
    print(pheatmap(mx,  color=color.heat3,scale="none",
         cluster_cols = F,
         breaks = seq(-2,2,length.out=101)))
dev.off()

