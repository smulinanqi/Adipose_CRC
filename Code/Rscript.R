colnames(rt)[2]<-"ftime"
colnames(rt)[1]<-"fustat"
my.surv <- Surv(rt$ftime, rt$fustat)

for (i in 3:ncol(rt)) {
  rt[,i]<-ifelse(rt[,i]>median(rt[,i]),'High','Low')
  
}

pl<-list()
for (i in 1:c(ncol(rt)-2)) {
  gene <-colnames(rt)[i+2]
  group = rt[,gene]
  data = cbind(rt[,1:2],group)
  group <- factor(group, levels = c('Low','High'))
  fit <- survfit(my.surv ~ group, data = data)
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  pl[[i]]<-ggsurvplot(fit, 
                      data = data,risk.table = T,
                      tables.theme = theme_cleantable(),
                      tables.height = 0.2,
                      fontsize=3,
                      tables.y.text=F,
                      legend=c(0.75,0.89),
                      legend.title =colnames(rt)[i+2],
                      legend.labs = c(paste0("High"," (",fit$n[1],")"),
                                      paste0("Low"," (",fit$n[2],")")),
                      xlab="Time (Month)",
                      ylab="DFS Rate for TCGA-COAD",
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"),
                      pval.size=3,
                      palette = c("#4682B4","#6E8B3D"),
                      censor=T,
                      censor.shape="|",
                      censor.size=2.5,
                      ggtheme=theme(axis.text.x = element_text(face="bold", color="black", size=8),    #各个字体大小
                                    axis.text.y = element_text(face="bold",  color="black", size=8),
                                    axis.title.x = element_text(face="bold", color="black", size=8),
                                    axis.title.y = element_text(face="bold",color="black", size=8),
                                    legend.text= element_text(face="bold", color="black", size=8),
                                    legend.title = element_text(face="bold", color="black", size=8),
                                    axis.line=element_line(color = "black"),
                                    # panel.border = element_blank(),
                                    # panel.grid =element_blank(),
                                    panel.background = element_blank(),
                                    plot.title=element_text(face="bold", color="black",size=8)))
}

length(pl)

res <- arrange_ggsurvplots(pl, 
                           print = F,
                           ncol = 3, nrow = 3)


tcga_expr<-df1
tcga_ddr<-df2
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_ddr)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_ddr[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}


data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data<-data[data$gene=='BMI',]



mut<-df
for (j in 3:4) {
  df<-mut[,c(1:2,j)]
  cor1<-df[,-ncol(df)]
  pathway<-colnames(cor1)
  LN_IC50<-as.numeric(df[,ncol(df)])
  
  box <- lapply(pathway,function(i) {
    dat=data.frame(Pathway=as.numeric(cor1[,i]),
                   group=LN_IC50) 
    head(dat)
    
    ggscatter(dat, x = "Pathway", y = "group",
              color = "black", shape = 19, size = 1, # Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "red", fill = "#BDBDBD"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
              xlab = i,
              legend = "right",
              ylab=colnames(mut)[j],
              cor.coeff.args = list(method = "spearman",  label.sep = "\n")
    )
  })
  plot_grid(plotlist=box, ncol=length(cor1) )
  ggsave(filename = paste0("Scatter/",colnames(mut)[j],".pdf"),width = 5,height = 2.26)
}


test<-rt
myUtest<-list()
for (j in colnames(test)[1:4]) {
  myUtest[[j]]<-lapply(colnames(test)[5:c(ncol(test))],function(i){
    y1=wilcox.test(test[,i]~test[,j],data=test)
    gene=i
    Pval=y1$p.value
    setNames(c(j ,gene, Pval ),
             c('G',"Group", "Pval"))
    
  })
  myUtest[[j]]<-do.call(rbind,myUtest[[j]])
}
myUtest<-do.call(rbind,myUtest)
myUtest<- as.data.frame(myUtest,
                        stringsAsFactors=FALSE)

for(i in c(3:3)){
  myUtest[,i] <- as.numeric(as.vector(myUtest[,i]))
}

trend<-test

dgm<-list()
for (j in colnames(trend)[c(5:ncol(trend))]) {
  dgm[[j]]<-as.data.frame(apply(as.data.frame(trend[,c(1:4)]), 2, function(x){
    tapply(trend[,j], INDEX = x, function(x) {mean(x)})##每一个纳入的病人对应的肿瘤亚型
  }))
}
dgm<-do.call(rbind,dgm)

ind<-unique(do.call(rbind,strsplit(rownames(dgm),'\\.'))[,1])

dgm1<-matrix(dimnames = list(ind,colnames(dgm)),nrow = length(ind),ncol = 4)
dgm1<-as.data.frame(dgm1)
for (i in 1:length(ind)) {
  dgm1[i,]<-dgm[i,]-dgm[i+1,]
}

for (i in 1:4) {
  dgm1[,i]<-ifelse(dgm1[,i]>0,'High-UP','High-Down')  
}

colnames(dgm1)<-gsub("\\-",'\\_',colnames(dgm1))
dgm1$Group<-rownames(dgm1)
dgm2<- dgm1%>%
  tidyr::gather(key="G",value="Trend",Best_ADI:Median_MUS) %>%
  dplyr::select(Group,G,Trend)

dgm2$G<-gsub('\\_',"\\-",dgm2$G) 
rownames(myUtest)<-paste0(myUtest$G,".",myUtest$Group)
rownames(dgm2)<-paste0(dgm2$G,'.',dgm2$Group)
dgm2<-dgm2[rownames(myUtest),]
identical(rownames(myUtest),rownames(dgm2))
mydf<-dgm2
mydf$Pval<-myUtest$Pval
colnames(mydf)[3]<-c("High_vs_Low")

mydf$High_vs_Low<-ifelse(mydf$Pval<0.0001,paste0(mydf$High_vs_Low,"****"),
                         ifelse(mydf$Pval<0.001,paste0(mydf$High_vs_Low,"***"),
                                ifelse(mydf$Pval<0.01,paste0(mydf$High_vs_Low,"**"),
                                       ifelse(mydf$Pval<0.05,paste0(mydf$High_vs_Low,"*"),
                                              paste0(mydf$High_vs_Low,'')))))
colnames(mydf)[1:2]<-c("ssGSEA","Group")


mydata1<-trend
for (i in 1:4) {
  mydata<-cbind(mydata1[,5:ncol(mydata1)],vale=mydata1[,i])
  mydata<-na.omit(mydata)
  mydata<- mydata%>%
    tidyr::gather(key="A",value="B",B_Cells_Memory:T_Cells_Regulatory_Tregs) %>%
    dplyr::select(vale,A,B)
  mydata$B<-mydata$B*100
  mydata$A<-gsub('\\_','\\ ',mydata$A)
  ggboxplot(mydata, 
            x = "A",
            y = "B",
            size=0.5,
            fill = "vale",
            ylab = "CIBERSORT (%)",
            xlab = "",
            add = "boxplot",title ='')+ 
    stat_compare_means(aes(group=vale),label = "p.signif",label.y=floor(max(mydata$B))+1,size=3.5,label.x = 1.3)+ theme(  ##format  signif
      axis.text.x = element_text(face="bold",  color="black", size=6,angle = 30,hjust=1),
      axis.text.y = element_text(face="bold",  color="black", size=6),
      axis.title.x = element_text(face="bold", color="black", size=6),
      axis.title.y = element_text(face="bold",color="black", size=6),
      legend.text= element_text(face="bold", color="black", size=6),
      legend.title = element_text(face="bold", color="black", size=6))+
    theme(legend.position = c(0.9,1.08),
          legend.background = element_rect(fill = "transparent"),
          legend.key.height =unit(0.1,"inches"),legend.key.width = unit(0.1,"inches"),
          legend.text = element_text(size = 6))+
    scale_fill_manual(labels = c("High","Low"),
                      values = c("#548B54","#DEB887"),name=do.call(rbind,strsplit(colnames(trend)[i],'\\-'))[,2])
  ggsave(filename = paste0("Box/",colnames(trend)[i],'_CIBERSORT_P_boxplot.pdf'),width = 7.4,height = 3.2)
}


forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),#
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))),
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#CD9B9B", lines="#8B7D7B", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(-1,0,1,2,3,4,5,6,7),
           lwd.xaxis=2,
           xlab=expression("log"[2]~"HR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "12" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.75,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
