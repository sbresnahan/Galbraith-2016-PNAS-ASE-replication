filter.SNPs <- function(counts,metadata,lcf){
  # Remove rows with < lcf counts counts by block and cross
  LA1 <- metadata[metadata$block==1&metadata$cross=="A","sample.id"]
  LA2 <- metadata[metadata$block==2&metadata$cross=="A","sample.id"]
  LB1 <- metadata[metadata$block==1&metadata$cross=="B","sample.id"]
  LB2 <- metadata[metadata$block==2&metadata$cross=="B","sample.id"]
  counts <- counts[rowSums(counts[,names(counts)%in%LA1])>lcf,]
  counts <- counts[rowSums(counts[,names(counts)%in%LA2])>lcf,]
  counts <- counts[rowSums(counts[,names(counts)%in%LB1])>lcf,]
  counts <- counts[rowSums(counts[,names(counts)%in%LB2])>lcf,]
  # Remove rows with < 2 SNPs
  counts$gene <- as.character(map(strsplit(row.names(counts), split = ":"), 2))
  genelist <- unique(counts$gene)
  delete.rows <- list()
  for(i in 1:length(genelist)){
    tmp <- counts[counts$gene==genelist[i],]
    if(length(row.names(tmp))<2){delete.rows <- append(delete.rows,genelist[i])}
  }
  counts <- counts[!counts$gene%in%unlist(delete.rows),]
  counts$gene <- NULL
  # Return filtered counts
  return(counts)
}

twobinom<-function(r1,n1,r2,n2,alpha=.05){
  # Modified for efficiency: changed outer() command to Rfast::Outer()
  n1p<-n1+1
  n2p<-n2+1
  n1m<-n1-1
  n2m<-n2-1
  q <- r1/n1
  p <- r2/n2
  if(is.na(q)){q <- 0}
  if(is.na(p)){p <- 0}
  chk<-abs(q-p)
  x<-c(0:n1)/n1
  y<-c(0:n2)/n2
  phat<-(r1+r2)/(n1+n2)
  m1<-t(Outer(x,y,"-"))
  m2<-matrix(1,n1p,n2p)
  flag<-(abs(m1)>=chk)
  m3<-m2*flag
  rm(m1,m2,flag)
  xv<-c(1:n1)
  yv<-c(1:n2)
  xv1<-n1-xv+1
  yv1<-n2-yv+1
  dis1<-c(1,pbeta(phat,xv,xv1))
  dis2<-c(1,pbeta(phat,yv,yv1))
  pd1<-NA
  pd2<-NA
  for(i in 1:n1){pd1[i]<-dis1[i]-dis1[i+1]}
  for(i in 1:n2){pd2[i]<-dis2[i]-dis2[i+1]}
  pd1[n1p]<-phat^n1
  pd2[n2p]<-phat^n2
  m4<-t(Outer(pd1,pd2,"*"))
  test<-sum(m3*m4)
  rm(m3,m4)
  list(p.value=test,p1=q,p2=p,est.dif=q-p)
}

PSGE.SK <- function(counts,metadata,phenotype){
  pat.exp <- counts[,metadata[metadata$parent%in%c("D")&
                                metadata$phenotype==phenotype,"sample.id"]]
  mat.exp <- counts[,metadata[metadata$parent%in%c("Q")&
                                metadata$phenotype==phenotype,"sample.id"]]
  i.len=length(row.names(pat.exp))
  return.df <- data.frame(matrix(ncol=2,nrow=0))
  names(return.df) <- c("SNP_gene","p")
  for(i in 1:i.len){
    print(i)
    SNP_gene=row.names(pat.exp[i,])
    p1.s=sum(pat.exp[i,])
    p2.s=sum(mat.exp[i,])
    p.o=sum(p1.s,p2.s)
    test=twobinom(r1=p1.s,n1=p.o,r2=p2.s,n2=p.o)$p.value
    return.df <- rbind(return.df, data.frame(SNP_gene=SNP_gene,p=test))
  }
  return(return.df)
}

PSGE.GLIMMIX <- function(counts,metadata){
  counts$SNP_gene <- row.names(counts)
  counts$geneID <- as.character(unlist(map(strsplit(counts$SNP_gene, 
                                                    split = ":"), 2)))
  genelist <- unique(counts$geneID)
  i.len <- length(genelist)
  df.out <- data.frame(matrix(ncol=4,nrow=0))
  names(df.out) <- c("ID","parent.p","cross.p","parentXcross.p")
  for(i in 1:i.len){
    counts.sub <- counts[counts$geneID==genelist[i],]
    counts.sub$geneID <- NULL
    counts.sub <- gather(counts.sub, sample.id, count, 
                         names(counts.sub), -SNP_gene, factor_key=TRUE)
    counts.sub <- join(counts.sub, metadata, by = "sample.id")
    counts.sub$parent <- as.factor(str_sub(counts.sub$parent,-1,-1))
    counts.sub$SNP_gene <- as.factor(counts.sub$SNP_gene)
    counts.sub$lineage <- as.factor(counts.sub$lineage)
    counts.sub$individual <- as.factor(counts.sub$individual)
    testfail <- F
    test <- "null"
    if(length(unique(counts.sub$SNP_gene))==1){
      tryCatchLog(test <- lmer(count~parent+lineage+parent*lineage+
                                 (1|individual), data=counts.sub), 
                  error = function(e) {testfail <- T})
    }else{
      tryCatchLog(test <- lmer(count~parent+lineage+parent*lineage+(1|SNP_gene)+
                                 (1|individual), data=counts.sub), 
                  error = function(e) {testfail <- T})}
    if(class(test)=="character"){testfail <- T}
    if(testfail==F){
      test <- summary(test)
      parent.p.list <- test[["coefficients"]][2,5]
      cross.p.list <- test[["coefficients"]][3,5]
      parent.cross.p.list <- test[["coefficients"]][4,5]
    }else{
      parent.p.list <- 1
      cross.p.list <- 1
      parent.cross.p.list <- 1
    }
    df.out <- rbind(df.out, data.frame(ID=genelist[i],
                                       parent.p=parent.p.list,
                                       cross.p=cross.p.list,
                                       parentXcross.p=parent.cross.p.list))
  }
  return(df.out)
}

PSGE.analysis <- function(counts,phenotype,metadata,SK,GLIMMIX){
  # Split count matrices by cross and parent of origin for plotting
  p1.pat <- counts[,metadata[metadata$parent%in%c("D")&
                               metadata$lineage=="EHB"&
                               metadata$phenotype==phenotype,"sample.id"]]
  p1.mat <- counts[,metadata[metadata$parent%in%c("Q")&
                               metadata$lineage=="EHB"&
                               metadata$phenotype==phenotype,"sample.id"]]
  p2.pat <- counts[,metadata[metadata$parent%in%c("D")&
                               metadata$lineage=="AHB"&
                               metadata$phenotype==phenotype,"sample.id"]]
  p2.mat <- counts[,metadata[metadata$parent%in%c("Q")&
                               metadata$lineage=="AHB"&
                               metadata$phenotype==phenotype,"sample.id"]]
  # Set up a data.frame to plot %p1 and %p2 for each SNP
  p1.plot <- data.frame(rowSums(p1.pat)/(rowSums(p1.mat)+rowSums(p1.pat)))
  names(p1.plot) <- c("p1")
  p1.plot[is.nan(p1.plot$p1),"p1"] <- 0
  p2.plot <- data.frame(rowSums(p2.mat)/(rowSums(p2.mat)+rowSums(p2.pat)))
  names(p2.plot) <- c("p2")
  p2.plot[is.nan(p2.plot$p2),"p2"] <- 0
  plot <- cbind(p1.plot,p2.plot)
  # Join results of Storer-Kim tests
  plot <- plot[row.names(plot)%in%SK$SNP_gene,]
  plot$SK.p <- SK$p
  plot$SNP_gene <- row.names(plot)
  plot$gene <- as.character(map(strsplit(plot$SNP_gene, split = ":"), 2))
  # Prep GLIMMIX results for downstream analyses
  GLIMMIX.biased <- data.frame(gene=GLIMMIX$ID,parent.p=GLIMMIX$parent.p,
                               cross.p=GLIMMIX$cross.p,parentXcross.p=GLIMMIX$
                                 parentXcross.p)
  # Correct for multiple testing
  plot$SK.padj <- p.adjust(plot$SK.p,"BH")
  plot$bias <- "NA"
  GLIMMIX$parent.padj <- p.adjust(GLIMMIX$parent.p,"BH")
  GLIMMIX$cross.padj <- p.adjust(GLIMMIX$cross.p,"BH")
  GLIMMIX$parentXcross.padj <- p.adjust(GLIMMIX$parentXcross.p,"BH")
  GLIMMIX.biased <- GLIMMIX[GLIMMIX$parent.padj<0.05|GLIMMIX$cross.padj<0.05,1]
  GLIMMIX.biased <- setdiff(GLIMMIX.biased,GLIMMIX[GLIMMIX$
                                                     parentXcross.padj<0.05,1])
  # For each gene, check whether all SNPs are biased in the same direction
  for(i in 1:length(row.names(plot))){
    p <- plot[i,"SK.padj"]
    p1 <- plot[i,"p1"]
    p2 <- plot[i,"p2"]
    if(p<0.05&p1>0.6&p2<0.4){plot[i,"bias"] <- "pat"}
    if(p<0.05&p1<0.4&p2>0.6){plot[i,"bias"] <- "mat"}
    if(p<0.05&p1<0.4&p2<0.4){plot[i,"bias"] <- "EHB"}
    if(p<0.05&p1>0.6&p2>0.6){plot[i,"bias"] <- "AHB"}
  }
  biaslist <- data.frame(matrix(ncol=2,nrow=0))
  names(biaslist) <- c("gene","bias")
  genelist <- unique(plot$gene)
  for(i in 1:length(genelist)){
    tmp <- unique(plot[plot$gene==genelist[i],"bias"])
    if(length(tmp)>1){
      if(length(tmp)==2){
        if(any(tmp%in%"NA")){
          bias <- tmp[!tmp%in%"NA"]
        }else{bias <- "NA"}
      }else{
        bias <- "NA"
      }
    }else{bias <- tmp}
    biaslist <- rbind(biaslist,data.frame(gene=genelist[[i]], bias=bias))
  }
  plot <- plot %>% left_join(biaslist, by = c('gene' = 'gene')) 
  names(plot)[c(7:8)] <- c("xbias","bias")
  plot$bias.plot <- "NA"
  for(i in 1:length(row.names(plot))){
    p1 <- plot$p1[i]
    p2 <- plot$p2[i]
    bias <- plot$bias[i]
    if(!bias=="NA"){
      if(bias=="pat"){if(p1>0.6&p2<0.4){plot[i,"bias.plot"]<- "pat"}}
      if(bias=="mat"){if(p1<0.4&p2>0.6){plot[i,"bias.plot"] <- "mat"}}
      if(bias=="EHB"){if(p1<0.4&p2<0.4){plot[i,"bias.plot"] <- "EHB"}}
      if(bias=="AHB"){if(p1>0.6&p2>0.6){plot[i,"bias.plot"] <- "AHB"}}
    }
  }
  plot[!plot$gene%in%GLIMMIX.biased,"bias.plot"] <- "NA" 
  plot <- rbind(plot[plot$bias.plot%in%c("NA"),],
                plot[plot$bias.plot%in%c("mat", "AHB", "EHB", "pat"),])
  plot$bias.plot <- factor(plot$bias.plot,
                           levels = c("NA","mat", "AHB", "EHB", "pat"))
  return(plot)
}

PSGE.collapse.avgExp <- function(data.plot,data.counts,metadata,phenotype){
  data.counts$SNP_gene <- row.names(data.counts)
  data <- data.plot %>% 
    left_join(data.counts, by = c('SNP_gene' = 'SNP_gene')) 
  genelist <- unique(data$gene)
  p1.mean <- list()
  p2.mean <- list()
  biaslist <- list()
  for(i in 1:length(genelist)){
    tmp <- data[data$gene==genelist[i],]
    if(!any(tmp$bias=="NA") & length(tmp[!tmp$bias.plot=="NA","p1"])>0){
      tmp.sub <- tmp[!tmp$bias.plot=="NA",]
      # Split count matrices by cross and parent of origin for plotting
      p1.pat <- tmp.sub[,metadata[metadata$parent%in%c("D")&
                                    metadata$lineage=="EHB"&
                                    metadata$phenotype==phenotype,"sample.id"]]
      p1.mat <- tmp.sub[,metadata[metadata$parent%in%c("Q")&
                                    metadata$lineage=="EHB"&
                                    metadata$phenotype==phenotype,"sample.id"]]
      p2.pat <- tmp.sub[,metadata[metadata$parent%in%c("D")&
                                    metadata$lineage=="AHB"&
                                    metadata$phenotype==phenotype,"sample.id"]]
      p2.mat <- tmp.sub[,metadata[metadata$parent%in%c("Q")&
                                    metadata$lineage=="AHB"&
                                    metadata$phenotype==phenotype,"sample.id"]]
      p1.mean.x <- mean(sum(p1.pat)/(sum(p1.mat)+sum(p1.pat)))
      if(is.nan(p1.mean.x)){p1.mean.x <- 0}
      if(is.infinite(p1.mean.x)){p1.mean.x <- 1}
      p1.mean[i] <- p1.mean.x
      p2.mean.x <- mean(sum(p2.mat)/(sum(p2.mat)+sum(p2.pat)))
      if(is.nan(p2.mean.x)){p2.mean.x <- 0}
      if(is.infinite(p2.mean.x)){p2.mean.x <- 1}
      p2.mean[i] <- p2.mean.x
      biaslist[i] <- as.character(tmp.sub$bias.plot[1])
    }else{
      p1.mean[i] <- mean(tmp$p1)
      p2.mean[i] <- mean(tmp$p2)
      biaslist[i] <- "NA"}}
  return.data <- data.frame(gene=unlist(genelist),
                            bias.plot=unlist(biaslist),
                            p1=unlist(p1.mean),
                            p2=unlist(p2.mean))
  return.data$bias.plot <- factor(return.data$bias.plot,
                                  levels = c("NA","mat", "AHB", "EHB", "pat"))
  return.data <- return.data[order(return.data$bias.plot),]
  return(return.data)
}

PSGE.plot.tx <- function(data,title){
  g <- ggplot(data, aes(x=p1, y=p2,
                        color=bias.plot,alpha=0.5)) + 
    geom_point(size=3) + theme_classic() +
    xlab(expression(paste("% A allele in ",E[mother],
                          " x ",A[father],sep=""))) +
    ylab(expression(paste("% A allele in ",A[mother],
                          " x ",E[father],sep=""))) +
    ggtitle(title) +
    theme(text = element_text(size=18),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(breaks = levels(data$bias.plot),
                       values=c("grey90","#2da330","#897a79","#c59638",
                                        "#e83a37")) +
                                          guides(alpha=F, color=F) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2))
  return(g)
}

triplot.plot <- function(sterile.plot,reproductive.plot,
                         label.A,label.B,allgenes){
  gmid.df <- data.frame(sterile=c(length(
    unique(sterile.plot[sterile.plot$bias.plot=="mat","gene"])),
    length(unique(sterile.plot[sterile.plot$bias.plot=="AHB","gene"])),
    length(unique(sterile.plot[sterile.plot$bias.plot=="EHB","gene"])),
    length(unique(sterile.plot[sterile.plot$bias.plot=="pat","gene"]))),
    Bias=c("mat","AHB","EHB","pat"),
    reproductive=c(length(unique(
      reproductive.plot[reproductive.plot$
                          bias.plot=="mat","gene"])),
      length(unique(reproductive.plot[reproductive.plot$bias.plot=="AHB","gene"])),
      length(unique(reproductive.plot[reproductive.plot$bias.plot=="EHB","gene"])),
      length(unique(reproductive.plot[reproductive.plot$bias.plot=="pat","gene"]))))
  
  ## Test if # of sterile biased genes is different from reproductive biased genes
  mat.test <- chisq.test(data.frame(Success=c(gmid.df[1,1],gmid.df[1,3]),
                                    Failure=c(length(allgenes)-gmid.df[1,1],
                                              length(allgenes)-gmid.df[1,3]),
                                    row.names=c("sterile",
                                                "reproductive")))$p.value
  AHB.test <- chisq.test(data.frame(Success=c(gmid.df[2,1],gmid.df[2,3]),
                                    Failure=c(length(allgenes)-gmid.df[2,1],
                                              length(allgenes)-gmid.df[2,3]),
                                    row.names=c("sterile",
                                                "reproductive")))$p.value
  EHB.test <- chisq.test(data.frame(Success=c(gmid.df[3,1],gmid.df[3,3]),
                                    Failure=c(length(allgenes)-gmid.df[3,1],
                                              length(allgenes)-gmid.df[3,3]),
                                    row.names=c("sterile",
                                                "reproductive")))$p.value
  pat.test <- chisq.test(data.frame(Success=c(gmid.df[4,1],gmid.df[4,3]),
                                    Failure=c(length(allgenes)-gmid.df[4,1],
                                              length(allgenes)-gmid.df[4,3]),
                                    row.names=c("sterile",
                                                "reproductive")))$p.value
  ## Build table
  gmid.df$`.` <- c(mat.test,AHB.test,EHB.test,pat.test)
  gmid.df <- gmid.df[,c(4,1,2,3)]
  nsrows <- row.names(gmid.df[gmid.df$`.`>0.05,])
  gmid.df$`.` <- formatC(gmid.df$`.`, format = "e", digits = 2)
  gmid.df[nsrows,"."] <- "(ns)"
  gmid.df <- gmid.df[,c(2,3,4,1)]
  cols <- matrix("black", nrow(gmid.df), ncol(gmid.df))
  cols[1,2] <- "#2da330"
  cols[2,2] <- "#897a79"
  cols[3,2] <- "#c59638"
  cols[4,2] <- "#e83a37"
  ccols <- matrix("white", nrow(gmid.df), ncol(gmid.df))
  ccols[1,3] <- "#f4efea"
  ccols[2,3] <- "#f4efea"
  ccols[3,3] <- "#f4efea"
  ccols[4,3] <- "#f4efea"
  ccols[1,1] <- "#f4efea"
  ccols[2,1] <- "#f4efea"
  ccols[3,1] <- "#f4efea"
  ccols[4,1] <- "#f4efea"
  ccols[1,2] <- "#e4d8d1"
  ccols[2,2] <- "#e4d8d1"
  ccols[3,2] <- "#e4d8d1"
  ccols[4,2] <- "#e4d8d1"
  cfonts <- matrix("plain", nrow(gmid.df), ncol(gmid.df))
  cfonts[1,2] <- "bold"
  cfonts[2,2] <- "bold"
  cfonts[3,2] <- "bold"
  cfonts[4,2] <- "bold"
  names(gmid.df) <- c(label.A,"Bias",label.B,".")
  gmid.df[2,2] <- "AHB"
  gmid.df[3,2] <- "EHB"
  tt <- ttheme_default(core=list(fg_params = list(col = cols, 
  cex = 1,
  fontface = cfonts),
  bg_params = list(col=NA,
  fill = ccols),
  padding.h=unit(2, "mm")),
  rowhead=list(bg_params = list(col=NA)),
  colhead=list(bg_params = list(fill = c("#f4efea",
  "#e4d8d1",
  "#f4efea",
  "white")),
  fg_params = list(rot=90,
  cex = 1,
  col = c("black",
  "black",
  "black",
  "white"))))
  
  gmid <- tableGrob(gmid.df, rows = NULL, theme=tt)
  return(gmid)
}