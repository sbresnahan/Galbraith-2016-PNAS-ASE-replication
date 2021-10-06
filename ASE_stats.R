library(purrr)
library(tidyverse)
library(viridis)
library(gridExtra)

SNPs <- read.table("SNPs_for_analysis.bed",sep="\t",header=F)[,c(1,2,4)]
names(SNPs) <- c("chr","pos","SNP_gene")
SNPs$SNP <- as.character(map(strsplit(SNPs$SNP_gene, split = ":"), 1))
SNPs$geneID <- as.character(map(strsplit(SNPs$SNP_gene, split = ":"), 2))
write.table(SNPs,"SNP_gene_overlaps.txt",sep="\t",quote=F,row.names=F)

SNP_counts <- data.frame(matrix(ncol=0,nrow=length(SNPs$SNP_gene)))
SNP_counts$SNP_gene <- SNPs$SNP_gene

files.counts <- list.files(path="counts_SNPs", 
                           pattern="*.txt", 
                           full.names=TRUE, 
                           recursive=FALSE)

for(i in 1:length(files.counts)){
  tmp <- read.table(files.counts[[i]],header=F)[,c(4,7)]
  tmp.name <-  strsplit(strsplit(files.counts[[i]], split = "/")[[1]][2],
                        split = "[.]")[[1]][1]
  tmp.name
  names(tmp) <- c("SNP_gene",tmp.name)
  SNP_counts <- SNP_counts %>% 
    left_join(tmp, by = c('SNP_gene' = 'SNP_gene')) 
}

row.names(SNP_counts) <- SNP_counts$SNP_gene
SNP_counts$SNP_gene <- NULL

SNP_counts[is.na(SNP_counts)] <- 0

write.csv(SNP_counts,"SNP_gene_counts.csv")

metadata <- read.csv("metadata.csv")

sterile_counts <- SNP_counts[,names(SNP_counts)%in%metadata[metadata$phenotype=="sterile","sample.id"]]
reproductive_counts <- SNP_counts[,names(SNP_counts)%in%metadata[metadata$phenotype=="reproductive","sample.id"]]

lcf <- 5

sterile_counts <- sterile_counts[rowSums(sterile_counts[,metadata[metadata$parent%in%c("875D","875Q")&metadata$phenotype=="sterile",
                                                                  "sample.id"]])>lcf,]

sterile_counts <- sterile_counts[rowSums(sterile_counts[,metadata[metadata$parent%in%c("888D","888Q")&metadata$phenotype=="sterile",
                                                                  "sample.id"]])>lcf,]

sterile_counts <- sterile_counts[rowSums(sterile_counts[,metadata[metadata$parent%in%c("882D","882Q")&metadata$phenotype=="sterile",
                                                                  "sample.id"]])>lcf,]

sterile_counts <- sterile_counts[rowSums(sterile_counts[,metadata[metadata$parent%in%c("894D","894Q")&metadata$phenotype=="sterile",
                                                                  "sample.id"]])>lcf,]

reproductive_counts <- reproductive_counts[rowSums(reproductive_counts[,metadata[metadata$parent%in%c("875D","875Q")&metadata$phenotype=="reproductive",
                                                                                 "sample.id"]])>lcf,]

reproductive_counts <- reproductive_counts[rowSums(reproductive_counts[,metadata[metadata$parent%in%c("888D","888Q")&metadata$phenotype=="reproductive",
                                                                                 "sample.id"]])>lcf,]

reproductive_counts <- reproductive_counts[rowSums(reproductive_counts[,metadata[metadata$parent%in%c("882D","882Q")&metadata$phenotype=="reproductive",
                                                                                 "sample.id"]])>lcf,]

reproductive_counts <- reproductive_counts[rowSums(reproductive_counts[,metadata[metadata$parent%in%c("894D","894Q")&metadata$phenotype=="reproductive",
                                                                                 "sample.id"]])>lcf,]


twobinom<-function(r1=sum(elimna(x)),n1=length(x),
                   r2=sum(elimna(y)),n2=length(y),
                   x=NA,y=NA,alpha=.05){
  #
  # Test the hypothesis that two independent binomials have equal
  # probability of success using the Storer--Kim method.
  #
  # Function from the WRS2 package
  #
  # r1=number of successes in group 1
  # n1=number of observations in group 1
  #
  n1p<-n1+1
  n2p<-n2+1
  n1m<-n1-1
  n2m<-n2-1
  chk<-abs(r1/n1-r2/n2)
  x<-c(0:n1)/n1
  y<-c(0:n2)/n2
  phat<-(r1+r2)/(n1+n2)
  m1<-outer(x,y,"-")
  m2<-matrix(1,n1p,n2p)
  flag<-(abs(m1)>=chk)
  m3<-m2*flag
  b1<-1
  b2<-1
  xv<-c(1:n1)
  yv<-c(1:n2)
  xv1<-n1-xv+1
  yv1<-n2-yv+1
  dis1<-c(1,pbeta(phat,xv,xv1))
  dis2<-c(1,pbeta(phat,yv,yv1))
  pd1<-NA
  pd2<-NA
  for(i in 1:n1)pd1[i]<-dis1[i]-dis1[i+1]
  for(i in 1:n2)pd2[i]<-dis2[i]-dis2[i+1]
  pd1[n1p]<-phat^n1
  pd2[n2p]<-phat^n2
  m4<-outer(pd1,pd2,"*")
  test<-sum(m3*m4)
  list(p.value=test,p1=r1/n1,p2=r2/n2,est.dif=r1/n1-r2/n2)
}

# Calculate p1 and p2 for each SNP
# Test each SNP for strength and direction of p1 and p2
test.df <- function(p1.pat.df,p1.mat.df,
                    p2.pat.df,p2.mat.df,high,low,
                    stest=c("CS","SK","thresholds")){
  return.list <- list()
  i.len <- length(row.names(p1.pat.df))
  for(i in 1:i.len){
    p1.pat <- p1.pat.df[i,]
    p1.mat <- p1.mat.df[i,]
    p2.pat <- p2.pat.df[i,]
    p2.mat <- p2.mat.df[i,]
    
    # p1 = sum(A in ExA)/(A in ExA + E in ExA)
    # p2 = sum(A in AxE)/(A in AxE + E in AxE)
    p1.s <- sum(p1.pat)
    p1.o <- sum(p1.pat,p1.mat)
    p1 <- p1.s/p1.o
    if(is.na(p1)){p1 <- 0}
    p2.s <- sum(p2.mat)
    p2.o <- sum(p2.pat,p2.mat)
    p2 <- p2.s/p2.o
    if(is.na(p2)){p2 <- 0}
    
    if(stest=="SK"){
      test <- twobinom(r1=p1.s,n1=p1.o,r2=p2.s,n2=p2.o)$p.value
    }
    
    if(stest=="CS"){
      test <- chisq.test(data.frame(X=c(p1.s,p2.s),Y=c(p1.o,p2.o)))$p.value
    }
    
    if(stest=="thresholds"){
      test <- 0
    }
    
    if(is.nan(test)){test <- 1}
    if(p1>high&p2<low){pat.test <- test}else{pat.test <- 1}
    if(p1<low&p2>high){mat.test <- test}else{mat.test <- 1}
    if(p1<low&p2<low){EHB.test <- test}else{EHB.test <- 1}
    if(p1>high&p2>high){AHB.test <- test}else{AHB.test <- 1}
    
    return.df <- data.frame(SNP_gene=row.names(p1.pat.df[i,]),
                            pat.p=pat.test,mat.p=mat.test,
                            AHB.p=AHB.test,EHB.p=EHB.test)
    
    return.list[[i]] <- return.df
  }
  return(bind_rows(return.list))
}

GLIMMIX <- function(counts,metadata){
  parent.list <- list()
  parent.p.list <- list()
  lineage.list <- list()
  lineage.p.list <- list()
  genelist <- unique(counts$ID)
  for(i in 1:length(genelist)){
    counts.sub <- counts[counts$ID==genelist[i],-c(1,3)]
    counts.sub <- gather(counts.sub, sample.id, count, 
                         names(counts.sub), -SNP.pos, factor_key=TRUE)
    counts.sub <- join(counts.sub, metadata, by = "sample.id")[,-c(4:6)]
    counts.sub$SNP.pos <- as.factor(counts.sub$SNP.pos)
    counts.sub$lineage <- as.factor(counts.sub$lineage)
    test <- summary(lme(count~parent.dir+lineage+parent.dir*lineage,
                        random=list(SNP.pos=~1),
                        data=counts.sub))
    parent.list[i] <- test$tTable[2,"Value"]
    parent.p.list[i] <- test$tTable[2,"p-value"]
    lineage.list[i] <- test$tTable[3,"Value"]
    lineage.p.list[i] <- test$tTable[3,"p-value"]
  }
  return(data.frame(ID=genelist,
                    parent.dirpat.est=unlist(parent.list),
                    parent.dirpat.p=unlist(parent.p.list),
                    lineageEHB.est=unlist(lineage.list),
                    lineageEHB.p=unlist(lineage.p.list)))
}



# Subset sterile worker SNP count matrix for input to sk.test.df function
## p1
p1.sterile.pat <- sterile_counts[,metadata[metadata$parent%in%c("875D","882D")&
                                             metadata$phenotype=="sterile","sample.id"]]
p1.sterile.mat <- sterile_counts[,metadata[metadata$parent%in%c("875Q","882Q")&
                                             metadata$phenotype=="sterile","sample.id"]]
names(p1.sterile.pat) <- names(p1.sterile.mat)
## p2
p2.sterile.pat <- sterile_counts[,metadata[metadata$parent%in%c("888D","894D")&
                                             metadata$phenotype=="sterile","sample.id"]]
p2.sterile.mat <- sterile_counts[,metadata[metadata$parent%in%c("888Q","894Q")&
                                             metadata$phenotype=="sterile","sample.id"]]
names(p2.sterile.pat) <- names(p2.sterile.mat)


sterile.SK <- test.df(p1.sterile.pat,p1.sterile.mat,
                      p2.sterile.pat,p2.sterile.mat,
                      low=0.4,high=0.6,"CS")


#sterile.GLIMMIX <- GLIMMIX(sterile_counts,metadata)


# Set up data.frame to plot p1 and p2 for each SNP
## p1 = sum(A in ExA)/(A in ExA + E in ExA)
p1.sterile.plot <- data.frame(rowSums(p1.sterile.pat)/(rowSums(p1.sterile.mat)+
                                                         rowSums(p1.sterile.pat)))
names(p1.sterile.plot) <- c("p1")

## p2 = sum(A in AxE)/(A in AxE + E in AxE)
p2.sterile.plot <- data.frame(rowSums(p2.sterile.mat)/(rowSums(p2.sterile.mat)+
                                                         rowSums(p2.sterile.pat)))
names(p2.sterile.plot) <- c("p2")

## Make data.frame from p1, p2, and sk.test.df output
sterile.plot <- cbind(p1.sterile.plot,p2.sterile.plot)
sterile.plot <- sterile.plot[row.names(sterile.plot)%in%sterile.SK$SNP_gene,]
sterile.plot <- cbind(sterile.plot,sterile.SK[,c(2:5)])
sterile.plot$SNP_gene <- row.names(sterile.plot)
sterile.plot$gene <- as.character(map(strsplit(sterile.plot$SNP_gene, split = ":"), 1))

## For each gene, check whether all SNPs are biased in the same direction
sterile.plot$bias <- "NA"
sterile.plot[sterile.plot$pat.p<0.05,"bias"] <- "pat"
sterile.plot[sterile.plot$mat.p<0.05,"bias"] <- "mat"
sterile.plot[sterile.plot$EHB.p<0.05,"bias"] <- "EHB"
sterile.plot[sterile.plot$AHB.p<0.05,"bias"] <- "AHB"
biaslist <- data.frame(matrix(ncol=2,nrow=0))
names(biaslist) <- c("gene","bias")
genelist <- unique(sterile.plot$gene)
for(i in 1:length(genelist)){
  tmp <- unique(sterile.plot[sterile.plot$gene==genelist[i],"bias"])
  if(length(tmp)>1){
    if(length(tmp)==2){
      if(any(tmp%in%"NA")){
        bias <- tmp[!tmp%in%"NA"]
      }
    }else{
      bias <- "NA"
    }
  }else{bias <- tmp}
  biaslist <- rbind(biaslist,data.frame(gene=genelist[[i]], bias=bias))
}
sterile.plot <- sterile.plot %>% 
  left_join(biaslist, by = c('gene' = 'gene')) 
names(sterile.plot) <- c("p1","p2","pat.p","mat.p",
                         "AHB.p","EHB.p","SNP_gene","gene","xbias","bias")

bias.plot <- list()
for(i in 1:length(row.names(sterile.plot))){
  if(sterile.plot[i,"xbias"]==sterile.plot[i,"bias"]){
    bias.plot[i] <- sterile.plot[i,"xbias"]
  }else{
    bias.plot[i] <- "NA"
  }
}
sterile.plot$bias.plot <- unlist(bias.plot)
sterile.plot <- rbind(sterile.plot[sterile.plot$bias.plot%in%c("NA"),],
                      sterile.plot[sterile.plot$bias.plot%in%c("mat", "AHB", "EHB", "pat"),])
sterile.plot$bias.plot <- factor(sterile.plot$bias.plot,
                                 levels = c("NA","mat", "AHB", "EHB", "pat"))

# Plot p1 and p2 for sterile workers
# rows (SNPs) of sterile.plot are colored in if the cooresponding gene is biased
g1 <- ggplot(sterile.plot, aes(x=p1, y=p2,
                               color=bias.plot,alpha=0.8)) + 
  geom_point() + theme_classic() +
  xlab(expression(paste("% A allele in ",E[mother],
                        " x ",A[father],sep=""))) +
  ylab(expression(paste("% A allele in ",A[mother],
                        " x ",E[father],sep=""))) +
  ggtitle("sterile workers") +
  theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = levels(sterile.plot$bias)[-c(1)],
                     values=c("grey", magma(20)[c(2)], 
                              magma(20)[c(8)], magma(20)[c(12)], 
                              magma(20)[c(16)]))+
  guides(alpha=F, color=F)


# Subset reproductive worker SNP count matrix for input to sk.test.df function
## p1
p1.reproductive.pat <- reproductive_counts[,metadata[metadata$parent%in%c("875D","882D")&
                                                       metadata$phenotype=="reproductive","sample.id"]]
p1.reproductive.mat <- reproductive_counts[,metadata[metadata$parent%in%c("875Q","882Q")&
                                                       metadata$phenotype=="reproductive","sample.id"]]
names(p1.reproductive.pat) <- names(p1.reproductive.mat)
## p2
p2.reproductive.pat <- reproductive_counts[,metadata[metadata$parent%in%c("888D","894D")&
                                                       metadata$phenotype=="reproductive","sample.id"]]
p2.reproductive.mat <- reproductive_counts[,metadata[metadata$parent%in%c("888Q","894Q")&
                                                       metadata$phenotype=="reproductive","sample.id"]]
names(p2.reproductive.pat) <- names(p2.reproductive.mat)


reproductive.SK <- test.df(p1.reproductive.pat,p1.reproductive.mat,
                           p2.reproductive.pat,p2.reproductive.mat,
                           low=0.4,high=0.6,"CS")


#reproductive.GLIMMIX <- GLIMMIX(reproductive_counts,metadata)


# Set up data.frame to plot p1 and p2 for each SNP
## p1 = sum(A in ExA)/(A in ExA + E in ExA)
p1.reproductive.plot <- data.frame(rowSums(p1.reproductive.pat)/(rowSums(p1.reproductive.mat)+
                                                                   rowSums(p1.reproductive.pat)))
names(p1.reproductive.plot) <- c("p1")

## p2 = sum(A in AxE)/(A in AxE + E in AxE)
p2.reproductive.plot <- data.frame(rowSums(p2.reproductive.mat)/(rowSums(p2.reproductive.mat)+
                                                                   rowSums(p2.reproductive.pat)))
names(p2.reproductive.plot) <- c("p2")

## Make data.frame from p1, p2, and sk.test.df output
reproductive.plot <- cbind(p1.reproductive.plot,p2.reproductive.plot)
reproductive.plot <- reproductive.plot[row.names(reproductive.plot)%in%reproductive.SK$SNP_gene,]
reproductive.plot <- cbind(reproductive.plot,reproductive.SK[,c(2:5)])
reproductive.plot$SNP_gene <- row.names(reproductive.plot)
reproductive.plot$gene <- as.character(map(strsplit(reproductive.plot$SNP_gene, split = ":"), 1))

## For each gene, check whether all SNPs are biased in the same direction
reproductive.plot$bias <- "NA"
reproductive.plot[reproductive.plot$pat.p<0.05,"bias"] <- "pat"
reproductive.plot[reproductive.plot$mat.p<0.05,"bias"] <- "mat"
reproductive.plot[reproductive.plot$EHB.p<0.05,"bias"] <- "EHB"
reproductive.plot[reproductive.plot$AHB.p<0.05,"bias"] <- "AHB"
biaslist <- data.frame(matrix(ncol=2,nrow=0))
names(biaslist) <- c("gene","bias")
genelist <- unique(reproductive.plot$gene)
for(i in 1:length(genelist)){
  tmp <- unique(reproductive.plot[reproductive.plot$gene==genelist[i],"bias"])
  if(length(tmp)>1){
    if(length(tmp)==2){
      if(any(tmp%in%"NA")){
        bias <- tmp[!tmp%in%"NA"]
      }
    }else{
      bias <- "NA"
    }
  }else{bias <- tmp}
  biaslist <- rbind(biaslist,data.frame(gene=genelist[[i]], bias=bias))
}
reproductive.plot <- reproductive.plot %>% 
  left_join(biaslist, by = c('gene' = 'gene')) 
names(reproductive.plot) <- c("p1","p2","pat.p","mat.p",
                              "AHB.p","EHB.p","SNP_gene","gene","xbias","bias")

bias.plot <- list()
for(i in 1:length(row.names(reproductive.plot))){
  if(reproductive.plot[i,"xbias"]==reproductive.plot[i,"bias"]){
    bias.plot[i] <- reproductive.plot[i,"xbias"]
  }else{
    bias.plot[i] <- "NA"
  }
}
reproductive.plot$bias.plot <- unlist(bias.plot)
reproductive.plot <- rbind(reproductive.plot[reproductive.plot$bias.plot%in%c("NA"),],
                           reproductive.plot[reproductive.plot$bias.plot%in%c("mat", "AHB", "EHB", "pat"),])
reproductive.plot$bias.plot <- factor(reproductive.plot$bias.plot,
                                      levels = c("NA","mat", "AHB", "EHB", "pat"))

g2 <- ggplot(reproductive.plot, aes(x=p1, y=p2,
                                    color=bias.plot,alpha=0.8)) + 
  geom_point() + theme_classic() +
  xlab(expression(paste("% A allele in ",E[mother],
                        " x ",A[father],sep=""))) +
  ylab(expression(paste("% A allele in ",A[mother],
                        " x ",E[father],sep=""))) +
  ggtitle("reproductive workers") +
  theme(text = element_text(size=20),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = levels(reproductive.plot$bias)[-c(1)],
                     values=c("grey", magma(20)[c(2)], 
                              magma(20)[c(8)], magma(20)[c(12)], 
                              magma(20)[c(16)]))+
  guides(alpha=F, color=F)


## Make tri-plot with #s of biased genes in middle section
gmid.df <- data.frame(Sterile=c(length(unique(sterile.plot[sterile.plot$bias=="mat","gene"])),
                                length(unique(sterile.plot[sterile.plot$bias=="AHB","gene"])),
                                length(unique(sterile.plot[sterile.plot$bias=="EHB","gene"])),
                                length(unique(sterile.plot[sterile.plot$bias=="pat","gene"]))),
                      Bias=c("mat","AHB","EHB","pat"),
                      Reproductive=c(length(unique(reproductive.plot[reproductive.plot$bias=="mat","gene"])),
                                     length(unique(reproductive.plot[reproductive.plot$bias=="AHB","gene"])),
                                     length(unique(reproductive.plot[reproductive.plot$bias=="EHB","gene"])),
                                     length(unique(reproductive.plot[reproductive.plot$bias=="pat","gene"]))))

sterile_IDs <- unique(as.character(map(strsplit(row.names(sterile_counts), split = ":"), 2)))
reproductive_IDs <- unique(as.character(map(strsplit(row.names(reproductive_counts), split = ":"), 2)))

## Test if # of sterile biased genes is different from reproductive biased genes
mat.test <- chisq.test(data.frame(Success=c(gmid.df[1,1],gmid.df[1,3]),
                                  Failure=c(length(sterile_IDs)-gmid.df[1,1],
                                            length(reproductive_IDs)-gmid.df[1,3]),
                                  row.names=c("Sterile","Reproductive")))$p.value
AHB.test <- chisq.test(data.frame(Success=c(gmid.df[2,1],gmid.df[2,3]),
                                  Failure=c(length(sterile_IDs)-gmid.df[2,1],
                                            length(reproductive_IDs)-gmid.df[2,3]),
                                  row.names=c("Sterile","Reproductive")))$p.value
EHB.test <- chisq.test(data.frame(Success=c(gmid.df[3,1],gmid.df[3,3]),
                                  Failure=c(length(sterile_IDs)-gmid.df[3,1],
                                            length(reproductive_IDs)-gmid.df[3,3]),
                                  row.names=c("Sterile","Reproductive")))$p.value
pat.test <- chisq.test(data.frame(Success=c(gmid.df[4,1],gmid.df[4,3]),
                                  Failure=c(length(sterile_IDs)-gmid.df[4,1],
                                            length(reproductive_IDs)-gmid.df[4,3]),
                                  row.names=c("Sterile","Reproductive")))$p.value
## Build table
gmid.df$`.` <- c(mat.test,AHB.test,EHB.test,pat.test)

gmid.df <- gmid.df[,c(4,1,2,3)]
nsrows <- row.names(gmid.df[gmid.df$`.`>0.05,])
gmid.df$`.` <- formatC(gmid.df$`.`, format = "e", digits = 2)
gmid.df[nsrows,"."] <- "(ns)"

cols <- matrix("black", nrow(gmid.df), ncol(gmid.df))
cols[1,3] <- magma(20)[c(2)]
cols[2,3] <- magma(20)[c(8)]
cols[3,3] <- magma(20)[c(12)]
cols[4,3] <- magma(20)[c(16)]

ccols <- matrix("white", nrow(gmid.df), ncol(gmid.df))
ccols[1,3] <- "#e4d8d1"
ccols[2,3] <- "#e4d8d1"
ccols[3,3] <- "#e4d8d1"
ccols[4,3] <- "#e4d8d1"
ccols[1,2] <- "#f4efea"
ccols[2,2] <- "#f4efea"
ccols[3,2] <- "#f4efea"
ccols[4,2] <- "#f4efea"
ccols[1,4] <- "#f4efea"
ccols[2,4] <- "#f4efea"
ccols[3,4] <- "#f4efea"
ccols[4,4] <- "#f4efea"

cfonts <- matrix("plain", nrow(gmid.df), ncol(gmid.df))
cfonts[1,3] <- "bold"
cfonts[2,3] <- "bold"
cfonts[3,3] <- "bold"
cfonts[4,3] <- "bold"

tt <- ttheme_default(core=list(fg_params = list(col = cols, 
                                                cex = 1,
                                                fontface = cfonts),
                               bg_params = list(col=NA,
                                                fill = ccols),
                               padding.h=unit(2, "mm")),
                     rowhead=list(bg_params = list(col=NA)),
                     colhead=list(bg_params = list(fill = c("white",
                                                            "#f4efea",
                                                            "#e4d8d1",
                                                            "#f4efea")),
                                  fg_params = list(rot=90,
                                                   cex = 1,
                                                   col = c("white",
                                                           "black",
                                                           "black",
                                                           "black"))))

gmid <- tableGrob(gmid.df, rows = NULL, theme=tt)

## Plot and save
fig1 <- grid.arrange(g1, gmid, g2, widths=c(5,2.5,5))
ggsave(file="fig1.pdf", plot=fig1, width=18, height=6)




