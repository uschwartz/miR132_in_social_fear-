setwd("/Volumes/PromisePegasus/Projects/12_20190626_miRNA_Array/analysis/result/")
res<-read.delim("conventional/res.IPscr_IPinh_extended.txt")
load("../obj/palmieri_final.rda")
library(stringr)
library(pheatmap)
library(Biobase)
library(RColorBrewer)
#positive in all 4
wx.target<-which(apply(res[,c("miRDB","miRMap","miRWalk2","targetScan")],1,sum)==4)

res.filt<-res[wx.target,c("SYMBOL","logFC","P.Value")]
res.filt[order(res.filt$P.Value),]

genes2cluster<-res[wx.target, 
                    grep(".CEL", colnames(res))]

rownames(genes2cluster)<-res[wx.target, "SYMBOL"]

genes2cluster


condition_names <- ifelse(str_detect(pData(palmieri_final)$condition,
                                     "input"), "input", "IP")
treatment_names <- ifelse(str_detect(pData(palmieri_final)$treatment,
                                     "scr"), "scr", "inh")

annotation_for_heatmap <-data.frame(condition = condition_names,  treatment = treatment_names)
row.names(annotation_for_heatmap) <- row.names(pData(palmieri_final))

anno_colors <- list( condition= c(input = "chartreuse4", IP = "burlywood3"),
                     treatment= c(inh = "blue4", scr = "cadetblue2"))


pheatmap(genes2cluster, annotation_col = annotation_for_heatmap,
         annotation_colors = anno_colors, scale = "row", show_rownames = F)


#####
cnames<-colnames(genes2cluster)
input<-apply(genes2cluster[,grep("_N_",cnames,value = T)],1,mean)
ip.scr<-apply(genes2cluster[,grep("_6F4_Scr",cnames,value = T)],1,mean)
ip.inh<-apply(genes2cluster[,grep("_6F4_Inh",cnames,value = T)],1,mean)

mx.pre<-cbind(input,ip.scr,ip.inh)

#get log(ratio)
row.mean<-apply(mx.pre,1,mean)
mx<-log(mx.pre/row.mean)

lv.samp<-colnames(mx)


# additional constrain up-regulated in IPs
# get enrichment in IP
mx.flt<-mx[(mx[,"input"]<mx[,"ip.scr"]),]


color.heat3 = colorRampPalette(rev(brewer.pal(n = 7,
                            name ="RdBu"))[c(rep(1,3),1:7,rep(7,3))])(100)
pdf("conventional/GDF5_heatmap.pdf", height=5, width = 4)
  pheatmap(mx.flt, cluster_cols=F,scale="none",
         color=color.heat3,
         breaks = seq(-0.15,0.15,length.out=101))
dev.off()
