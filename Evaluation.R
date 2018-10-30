library(plyr)
library(dplyr)
library(reshape2)
library(mclust)
library(NMI)
library(clv)



MatrixEuclidean <- function(data){
  smat <- apply(data, 1, crossprod)
  mat1 <- matrix(smat, nrow=nrow(data),ncol = nrow(data))
  mat3 <- tcrossprod(data)
  mat4 <- mat1 + t(mat1) - 2*mat3
  diag(mat4) <- 0
  return(sqrt(mat4))
}

randIndex <-function(group1,group2){
  group1 = unlist(group1)
  group2 = unlist(group2)
  x <- sapply(group1, function(x) x == group1)
  y <- sapply(group2, function(x) x == group2)
  sg <- sum(!(x == y))/2
  bc <- choose(dim(x)[1], 2)
  ri <- 1 - sg/bc
  return(ri)
}

align <- function(data.matrix){
  data.matrix = as.matrix(data.matrix)
  aligned = data.frame()
  temp = data.matrix
  
  while(1){
    r = which(temp == max(temp),arr.ind = T)
    rname = row.names(temp)
    cname = colnames(temp)
    rn = data.frame(row = c(rname[r[1,1]]),col = c(cname[r[1,2]]))
    aligned = rbind(aligned,rn)
    temp = data.matrix[setdiff(rname,aligned$row),setdiff(cname,aligned$col)]
    if(is.matrix(temp)){}
    else{
      aligned = rbind(
        aligned,
        data.frame(
          row = setdiff(rname,aligned$row),
          col = setdiff(cname,aligned$col)
        ))
      break
    }
  }
  row.names(aligned) = aligned$col
  ordered = data.matrix[as.character(aligned[colnames(data.matrix),]$row),]
  return(list(aligned.lable = aligned,aligned.matrix = ordered))
}


count2Jaccard <- function(data){
  rs = rowSums(data)
  cs = colSums(data)
  denominator =  sapply(cs,function(x){rs+x})-data
  return(data/denominator)
}

count2Fscore <- function(data,label = 'col'){
  rs = rowSums(data)
  cs = colSums(data)
  c = diag(data)
  Fscore = 2.0/(rs/c+cs/c)
  if(label == 'col'){
    names(Fscore) = colnames(data)
  }
  else{
    names(Fscore) = rownames(data)
  }
  return(Fscore)
}


evaluateClusterings <- function(label,true.label=NULL,data = NULL,n.dims.to.show = 4){
  
  external.eval = list()
  internal.eval = list()
  mcols = c('#2BBDC2','#F27E22','#34B64B','#7191CC','red','#C96BAC','#F7D729','#3DAAE1','#602390','#D7006F')
  hcols = colorRampPalette(c('white','red'))(101)  
  label = unlist(label)
  if(is.null(true.label) && is.null(data)){
    print('Error: at least one of true.label and data should be provided ')
    return()
  }
  if(!is.null(true.label)){
    true.label = unlist(true.label)
    if(length(label)!=length(true.label)){
      print('Error: two label vectors should have the same length') 
      return()
    }
    ARI = adjustedRandIndex(label,true.label)
    RI = randIndex(label,true.label)
    NMI <- NMI(data.frame(1:length(label), label), data.frame(1:length(true.label), true.label))
    m = table(label,true.label)
    J.m = count2Jaccard(m)
    J.m = align(J.m)$aligned.matrix
    aligned.label = align(J.m)$aligned.lable
    rownames(aligned.label) = aligned.label$col
    F.m = count2Fscore(m[as.character(aligned.label[colnames(m),1]),])
    F.m = diag(F.m)
    m = m[as.character(aligned.label[colnames(m),1]),]
    row.names(J.m) = paste('C',1:6,sep = '')
    colnames(F.m) = names(table(label))
    row.names(F.m) = paste('C',1:6,sep = '')
    row.names(m) = paste('C',1:6,sep = '')
    J.m.df = melt(J.m)
    F.m.df = melt(F.m)
    m.df = melt(m)
    Jaccard.plot <- ggplot(J.m.df,aes(true.label,label[nrow(J.m.df):1])) + geom_tile(aes(fill = value)) + 
      scale_fill_continuous(
        low = hcols[round(min(J.m.df$value),digits = 2)*100+1],
        high = hcols[round(max(J.m.df$value),digits = 2)*100+1]) + 
      geom_text(aes(label = round(value,digits = 2)),size = 3) + labs(x = 'Originl clusters',y = 'Reproduced clusters') +
      theme_bw() + theme(
        legend.position = 'none',
        axis.title = element_text(size = 8,family = 'Helvetica',face = 'bold'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      )
    Fscore.plot <- ggplot(F.m.df,aes(Var2,Var1[nrow(J.m.df):1])) + geom_tile(aes(fill = value)) + 
      scale_fill_continuous(
        low = hcols[round(min(F.m.df$value),digits = 2)*100+1],
        high = hcols[round(max(F.m.df$value),digits = 2)*100+1]) + 
      geom_text(aes(label = round(value,digits = 2)),size = 3) + labs(x = 'Originl clusters',y = 'Reproduced clusters') +
      theme_bw() + theme(
        legend.position = 'none',
        axis.title = element_text(size = 8,family = 'Helvetica',face = 'bold'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      ) 
    
     m.plot <- ggplot(m.df,aes(true.label,label[nrow(m.df):1])) + geom_tile(aes(fill = value)) + 
      scale_fill_continuous(
        low = 'white',
        high = 'red') + 
      geom_text(aes(label = round(value,digits = 2)),size = 3) + labs(x = 'Originl clusters',y = 'Reproduced clusters') + 
      theme_bw() + theme(
        legend.position = 'none',
        axis.title = element_text(size = 8,family = 'Helvetica',face = 'bold'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
      ) 
    external.eval = list(
      aligned.label = aligned.label,
      ARI = ARI,
      RI = RI,
      NMI = NMI$value,
      overall.Fscore = mean(diag(F.m)),
      cluster.Fscore = diag(F.m),
      confusion.map = m,
      Jaccard.map = J.m,
      Jaccard.plot=Jaccard.plot,
      Fscore.plot = Fscore.plot,
      m.plot = m.plot
      
    )
  }
  if(!is.null(data)){
    label1 = mapvalues(label,from = levels(factor(label)),to = 1:length(unique(label)))
    cls <- cls.scatt.data(data, as.integer(as.vector(label1)), dist="euclidean")
    dunn_temp <- clv.Dunn(cls, "complete", "complete")
    dunn.index = dunn_temp[1,1]
    
    DIS <- MatrixEuclidean(data)
    SIL <- silhouette(as.integer(as.vector(label1)), dmatrix = DIS)
    cluster.silhouette = tapply(SIL[,3],label,mean)
    overall.silhouette = mean(cluster.silhouette)
    S = as.data.frame(SIL[,1:3])
    S = setorder(S,-cluster,sil_width)
    S$ID = 1:nrow(S)
    silhouette.plot <- ggplot(as.data.frame(S),aes(ID,sil_width)) + 
      geom_bar(aes(color = as.factor(cluster),
                   fill = as.factor(cluster)),
               stat = 'identity',
               position = position_dodge()) +  coord_flip() + ylim(min(S[,3]),1) + 
      theme_classic() + theme(
        legend.position = 'none',
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) + scale_fill_manual(values = mcols) + scale_color_manual(values = mcols)  
    #      geom_text(aes(x = rep(1,length(unique(SIL[,1]))), y = as.data.frame(tapply(S$ID,S$cluster,mean)),label = 1:6))
    
    df =cbind(
      as.data.frame(data[,1:n.dims.to.show]),
      label = factor(label,levels = as.character(aligned.label[order(as.character(aligned.label$col)),]$row))
    )
    p <- ggpairs(
      df,
      columns = 1:n.dims.to.show,
      lower = list(
        continuous = 'points',
        mapping = aes(colour = label,fill = label),
        size = 0.5),
      diag = list(mapping = aes(colour = label,fill = label,alpha = 0.2)),
      upper = list(
        continuous = 'points',
        mapping = aes(colour = label,fill = label ),
        size = 0.5)
    )  + theme_bw() + theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 20) 
    )      
    for(t in 1:p$nrow){
      for(s in 1:p$ncol){
        p[t,s] = p[t,s] + scale_color_manual(values = mcols) + scale_fill_manual(values = mcols)
      }
    }
    
    internal.eval = list(
      dunn.index = dunn.index,
      overall.silhouette = overall.silhouette,
      sample.silhouette = S,
      cluster.silhouette = cluster.silhouette,
      silhouette.plot = silhouette.plot,
      dim.plot = p
    )
  }
  
  return(list(external.eval = external.eval,internal.eval = internal.eval))
}






