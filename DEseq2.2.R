#!/usr/bin/Rscript
suppressPackageStartupMessages({
  library("data.table")
  library(DESeq2)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(ggpubr)
  library(ggthemes)
  library("tibble")
  library("optparse")
  library("BiocParallel")
  
})

option_list <- list(
  make_option(c("-i", "--inputfile"), dest = "Rinput", default = "", 
              help="[opt] inputfile with onle reads counts"),
  make_option(c("-n", "--namefle"), dest = "namefile", default = "",
              help="[opt] namefile"),
  make_option(c("-c","--condition"),dest = "conditonfile",default = "",
              help = "[opt] condition config file"),
  make_option(c("-t","--contrast"),dest = "contrastfile",default = "",
              help = "[opt] contrast config file"),
  make_option(c("-p", "--padj"), dest = "p", default = "",
              help="[opt] padj cutoff "),
  make_option(c("-f", "--foldchange"), dest = "f", default = "0",
              help="[opt] log2foldchange defalut 0"),
  make_option(c("-x", "--threads"), dest = "x", default = "4",
              help="[opt] threads defalut 4"),
  make_option(c("-o", "--outprefix"), dest = "outfile", default = "",
              help="[opt] DEseq2 out file")
) 

parser <- OptionParser(usage = "%prog [options] file",
     option_list=option_list, description = " 2019-3-24 \
Description: \
  Wrapper for DEseq2.\
Example: \
  ./DEseq2.R -i ~/project/TE_ann/featureCountOUT/split/AABBRinput.txt -n name.conf -c condition.conf -t contrast.conf -o ./testet -p 1 -f 0 \ "
)
#
arguments <- parse_args(parser)
opt <- arguments$options
#infile
test <- arguments$Rinput
#test <- "/data/user/qinz/software/DEseq2/YQW/YQW.marix"
####多核计算
register(MulticoreParam(arguments$x))

if(test == "") {
  print_help(parser)
  stop(sprintf("input file is not specified"))
  test = file("stdin")  
}

namefile <- arguments$namefile
#namefile <- "/data/user/qinz/software/DEseq2/YQW/name.conf"
#namefile <- "~/qinz/TXJ/190777-16-17/data/name.conf"
conditonfile <- arguments$conditonfile
#conditonfile <- "/data/user/qinz/software/DEseq2/YQW/condition.conf"
contrastfile <- arguments$contrastfile
#contrastfile <- "/data/user/qinz/software/DEseq2/YQW/contrast.conf"



P <- as.numeric(arguments$p)
#p = 0.05
F <- as.numeric(arguments$f)
#F = 1
o <- arguments$outfile
#o <- "YQW_DEseq"

#workdir<-getwd()

test2 = read.csv(test,sep = "\t",header = FALSE)
namefile = read.csv(namefile,sep = "\t",header = FALSE)
conditonfile = read.csv(conditonfile,sep="\t",header = TRUE)
contrastfile = read.csv(contrastfile,sep="\t",header = FALSE)

smatrix <- as.matrix(test2[2:length(test2)])

resampleNames <- as.vector(namefile$V1)
row.names(smatrix) <- test2[,1]
colnames(smatrix) <- resampleNames

colData <- data.frame(row.names = colnames(smatrix),conditonfile)
dds <-DESeqDataSetFromMatrix(smatrix,colData = colData,design= ~ condition)
#dds <- dds[ rowSums(counts(dds)) > 0, ]
dds2 <- DESeq(dds,parallel = TRUE)

#rld <- rlog(dds2)
#
mergetable <- data.frame(row.names =row.names(dds2))
 
for (i in 1:nrow(contrastfile)) {
    contrast_name <- paste0(as.character(unlist(contrastfile[i,]))[2:3],collapse = "vs")
    contrast_vector <- as.vector(unlist(contrastfile[i,]))
    results <- results(dds2, contrast = contrast_vector,parallel = T)
    results <- lfcShrink(dds2, contrast = contrast_vector, res=results,parallel = T)
    if( arguments$p  == "") {
        diff <- subset(results,(log2FoldChange > F|log2FoldChange < -F)) %>% as.data.frame()
    }else {
    diff <- subset(results, padj<P & (log2FoldChange > F|log2FoldChange < -F)) %>% as.data.frame()
    }
    diff2 <- subset(results) %>% as.data.frame()
    diff2$logQ <- -log10(diff2$padj)
    diff2$Group <- "Not-change"
    diff2$Label <- ""
    diff2 <- diff2[order(diff2$padj),] 
    diff2$Group[which( (diff2$padj<0.05) & (diff2$log2FoldChange > F))] = "Up-regulated"
    diff2$Group[which( (diff2$padj<0.05) & (diff2$log2FoldChange < -F) )] = "Down-regulated"
    upgenes <- head(rownames(diff2[which(diff2$Group =="Up-regulated"),]),10)
    downgenes <- head(rownames(diff2[which(diff2$Group =="Down-regulated"),]),10)
    top10genes <- c(as.character(upgenes),as.character(downgenes))
    diff2$Label[match(top10genes,rownames(diff2)) ] <- top10genes
    Vplot <- ggscatter(diff2,x = "log2FoldChange",y = "logQ", color = "Group" , palette=c("#2f5688","#BBBBBB","#CC0000"),size=1,label =diff2$Label,font.label = 8,repel = T,xlab = "Log2FC",ylab="-Log10(padj)") + theme_base() +
      geom_hline(yintercept = -log10(P), linetype="dashed")+
      geom_vline(xintercept = -c(-F, F), linetype="dashed")
    ggsave(paste0("./",contrast_name,"vp.pdf"),plot = Vplot)
    colnames(diff) <-  paste(contrast_name,colnames(diff),sep = "_")
    mergetable <- merge(mergetable,diff[,c(2,6)],by="row.names",sort=FALSE,all=TRUE)
    rownames(mergetable) <- mergetable$Row.names
    mergetable$Row.names <- NULL
    final <- merge(diff,as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    outfile <- paste(o,contrast_name,sep = ".")
    write.table(final,file = outfile,row.names = FALSE,sep = "\t",quote = FALSE)
}

mergetable <- merge(mergetable,counts(dds2,normalize=TRUE),by="row.names",sort=FALSE,all=TRUE)

write.table(mergetable,file = paste0(o,"all"),row.names = FALSE,sep = "\t",quote = FALSE)


