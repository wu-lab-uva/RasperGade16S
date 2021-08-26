#' @export
calculate_Aitchison_dist = function(abundance.table){
  if(class(abundance.table)!="phyloseq") abundance.table = phyloseq(otu_table(object = abundance.table,taxa_are_rows = FALSE))
  clr = microbiome::transform(x = abundance.table,transform = "clr")
  DM <- vegdist(x = clr, method = "euclidean")
  return(DM)
}
#' @export
calculateACN = function(GCN,abundance,as.gene.abundance=FALSE){
  if(as.gene.abundance){
    acn = sum(abundance)/sum(abundance/GCN)
  }else{
    acn = sum(abundance*GCN)/sum(abundance)
  }
  return(acn)
}
#' @export
calculateACN_CI_bySimulation = function(abundance,error,discrete=FALSE,n=1000,alpha=0.05,
                                        detail=FALSE,as.gene.abundance=FALSE){
  sim.ACN = sapply(1:n,function(nn){
    this.GCN = sapply(error,function(x){rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)})
    if(discrete) this.GCN = round(this.GCN)
    this.GCN[this.GCN<=1] = 1
    calculateACN(GCN = this.GCN,abundance = abundance,as.gene.abundance = as.gene.abundance)
  })
  if(detail) return(sim.ACN)
  return(quantile(sim.ACN,probs = c(alpha/2,1-alpha/2)))
}
#' @export
classify_oligotrphic_community_by_abundance = 
  function(abundance,GCN,GCN.threshold = 2.5,abundance.threshold = c(0.4,0.6),as.gene.abundance=FALSE){
    if(as.gene.abundance) abundance = abundance/GCN
    oligo.proportion = sum(abundance[GCN<GCN.threshold])/sum(abundance)
    return(sum(as.numeric(oligo.proportion<abundance.threshold))-1)
  }
#' @export
classify_oligotrphic_community_by_abundance_with_confidence = 
  function(abundance,error,n=1000,
           GCN.threshold = 2.5,abundance.threshold = c(0.4,0.6),as.gene.abundance=FALSE){
    sim.oligo = sapply(1:n,function(nn){
      this.GCN = sapply(error,function(x){rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)})
      this.GCN = round(this.GCN)
      this.GCN[this.GCN<=1] = 1
      classify_oligotrphic_community_by_abundance(
        abundance = abundance,GCN = this.GCN,GCN.threshold = GCN.threshold,
        abundance.threshold = abundance.threshold,as.gene.abundance = as.gene.abundance)
    })
    oligo.freq = c(oligo=sum(sim.oligo<0),meso=sum(sim.oligo==0),copio=sum(sim.oligo>0))
    max.freq = max(oligo.freq)
    type = which(oligo.freq==max.freq)-2
    if(length(type)>1) type = type[which.min(abs(type))]
    return(list(type = unname(type),confidence = unname(max.freq/n),detail = oligo.freq))
  }
#' @export
classify_oligotrphic_community_by_ACN =
  function(abundance,GCN,as.gene.abundance=FALSE,threshold=3){
    (sum(calculateACN(GCN = GCN,abundance = abundance,as.gene.abundance = as.gene.abundance)>threshold)-
       length(threshold)/2)*(2/length(threshold))
  }
#' @export
classify_oligotrphic_community_by_ACN_with_confidence =
  function(abundance,error,as.gene.abundance=FALSE,threshold=3,n=1000,detail=FALSE){
    sim.ACN = calculateACN_CI_bySimulation(abundance = abundance,error = error,n = n,detail = TRUE,
                                           as.gene.abundance = as.gene.abundance,discrete = TRUE)
    this.classify = sapply(sim.ACN,function(x){
      (sum(x>threshold)-length(threshold)/2)*(2/length(threshold))
    })
    oligo.freq = c(oligo=sum(this.classify<0),meso=sum(this.classify==0),copio=sum(this.classify>0))
    max.freq = max(oligo.freq)
    type = which(oligo.freq==max.freq)-2
    if(length(type)>1) type = type[which.min(abs(type))]
    if(detail) return(list(type = unname(type),confidence = unname(max.freq/n),
                           detail = oligo.freq,ACN=sim.ACN))
    return(list(type = unname(type),confidence = unname(max.freq/n),detail = oligo.freq))
  }
#' @export
correct_RAD = function(abundance,GCN){
  this.rad = sort(relative_abundance(abundance/round(GCN)),decreasing = TRUE)
  return(this.rad)
}
#' @export
correct_RAD_with_confidence = function(abundance,error,lower.edge="bounce",n=1000,detail=FALSE){
  perm.GCN = sapply(1:n,function(nn){
    this.GCN = sample_GCN_from_error_distribution(error = error,lower.edge = lower.edge)
    return(this.GCN)
  })
  perm.abun = sapply(1:n,function(nn){
    this.perm.abun = abundance/round(perm.GCN[,nn])
    return(this.perm.abun)
  })
  perm.rad = lapply(1:n,function(nn){
    this.abun = perm.abun[,nn]
    names(this.abun) = names(abundance)
    this.rad = sort(relative_abundance(this.abun),decreasing = TRUE)
    return(this.rad)
  })
  median.abun = apply(perm.abun,1,median)
  CI.abun = apply(apply(perm.abun,2,relative_abundance),1,quantile,probs=c(0.025,0.975))
  names(median.abun) = names(abundance)
  colnames(CI.abun) = names(abundance)
  this.rad = sort(relative_abundance(median.abun),decreasing = TRUE)
  rank.name = sapply(perm.rad,names)
  rank.support = sapply(1:length(this.rad),function(this.rank){
    local.support = sum(rank.name[this.rank,]==names(this.rad)[this.rank])/n
    accumulative.support = sum(unlist(as.vector(rank.name[1:this.rank,]))==names(this.rad)[this.rank])/n
    return(c(local=local.support,accumulate = accumulative.support))
  })
  raw.support = sapply(1:length(abundance),function(this.rank){
    this.abundance = sort(abundance,decreasing = TRUE)
    local.support = sum(rank.name[this.rank,]==names(this.abundance)[this.rank])/n
    accumulative.support = sum(unlist(as.vector(rank.name[1:this.rank,]))==names(this.abundance)[this.rank])/n
    return(c(local=local.support,accumulate = accumulative.support))
  })
  if(detail) return(list(rad=this.rad,support = rank.support,
                         CI=CI.abun[,names(this.rad)],gene.support=raw.support,
                         GCN = perm.GCN,abundance=perm.abun))
  return(list(rad=this.rad,support = rank.support,CI=CI.abun[,names(this.rad)],gene.support=raw.support))
}
#' @export
coverage_length = function(rad,CI,detail=FALSE){
  cvl = sapply(1:length(rad),function(i){
    (rad[i]<=CI[2,i])&(rad[i]>=CI[1,i])
  })
  if(detail) return(cvl)
  return(sum(cvl)/length(rad))
}
#' @export
rank_square_difference = function(rank1,rank2,average=TRUE,normalize=FALSE){
  N= length(rank1)
  sum.diff = sum((1:N-as.numeric(factor(x = rank1,levels = rank2)))^2)
  max.diff = N*(N+1)*(N-1)/3
  if(normalize){
    sum.diff = -2*sum.diff/max.diff+1
  }else{
    if(average) sum.diff = sum.diff/N
  }
  return(unname(sum.diff))
}
#' @export
mean_support = function(support){
  apply(support,1,function(x){exp(mean(log(x)))})
}
#' @export
square_difference = function(x,y,average=TRUE){
  if(average) return(mean((x-y)^2))
  return((x-y)^2)
}
#' @export
correct_distance_matrix = function(abundance.table,GCN,method="bray",...){
  GCN = round(GCN)
  correct.table = correct_abundance_table(abundance.table,GCN)
  correct.phyloseq = phyloseq(otu_table(correct.table,taxa_are_rows = FALSE))
  if(method=="Aitchison") return(calculate_Aitchison_dist(correct.phyloseq))
  return(distance(physeq = correct.phyloseq,method = method,...))
}
#' @export
correct_distance_matrix_with_confidence = function(abundance.table,error,method="bray",alpha=0.05,
                                                   dimension=dim(abundance.table)[1],n=1000,detail=FALSE,numCores=1,...){
  perm.GCN = sapply(1:n,function(nn){
    this.GCN = sapply(error,function(x){rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)})
    this.GCN[this.GCN<=1] = 1
    return(this.GCN)
  })
  GCN = round(apply(perm.GCN,1,median))
  correct.table = correct_abundance_table(abundance.table,GCN)
  correct.phyloseq = phyloseq(otu_table(correct.table,taxa_are_rows = FALSE))
  if(method=="Aitchison"){
    correct.dist = calculate_Aitchison_dist(correct.phyloseq)
  }else{
    correct.dist = distance(physeq = correct.phyloseq,method = method,...)
  }
  perm.dist = mclapply(1:dim(perm.GCN)[2],function(k){
    this.GCN = round(perm.GCN[,k])
    this.table = correct_abundance_table(abundance.table,this.GCN)
    this.phyloseq = phyloseq(otu_table(this.table,taxa_are_rows = FALSE))
    if(method=="Aitchison"){
      this.dist = calculate_Aitchison_dist(this.phyloseq)
    }else{
      this.dist = distance(physeq = this.phyloseq,method = method,...)
    }
    return(this.dist)
  },mc.cores = numCores)
  alpha = 1-exp(log(1-alpha)/(dimension-1))
  CI.dist = sapply(1:length(perm.dist[[1]]),function(j){
    quantile(sapply(perm.dist,function(x){x[j]}),probs = c(alpha/2,1-alpha/2))
  })
  if(detail) return(list(dist = correct.dist,CI=CI.dist,perm=perm.dist))
  return(list(dist = correct.dist,CI=CI.dist))
}
#' @export
total_shift = function(dist1,dist2,squared=TRUE){
  if(squared){
    sum(apply(as.matrix(dist1-dist2),1,function(x){sum(x^2)}))
  }else{
    sum(apply(as.matrix(dist1-dist2),1,function(x){sqrt(sum(x^2))}))
  }
}
#' @export
coverage_dist = function(dist1,CI,detail=FALSE){
  uCI = as.matrix(CI[2,]+dist1*0)
  lCI = as.matrix(CI[1,]+dist1*0)
  mat = as.matrix(dist1)
  coverage = sapply(1:dim(mat)[1],function(i){
    all((mat[i,-i]<=uCI[i,-i])&(mat[i,-i]>=lCI[i,-i]))
  })
  if(detail) return(coverage)
  return(sum(coverage)/length(coverage))
}
#' @export
elementwise_coverage_dist = function(dist1,CI,detail=FALSE){
  uCI = as.matrix(CI[2,]+dist1*0)
  lCI = as.matrix(CI[1,]+dist1*0)
  mat = as.matrix(dist1)
  coverage = (mat<=uCI)&(mat>=lCI)
  if(detail) return(coverage)
  return(sum(coverage)/prod(dim(coverage)))
}
#' @export
correct_PERMANOVA_with_confidence = function(abundance.table,error,metadata,
                                             method="bray",alpha=0.05,n=1000,numCores=1,...){
  perm.GCN = sapply(1:n,function(nn){
    this.GCN = sapply(error,function(x){rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)})
    this.GCN[this.GCN<=1] = 1
    return(this.GCN)
  })
  GCN = round(apply(perm.GCN,1,median))
  correct.table = correct_abundance_table(abundance.table,GCN)
  correct.phyloseq = phyloseq(otu_table(correct.table,taxa_are_rows = FALSE))
  if(method=="Aitchison"){
    correct.dist = calculate_Aitchison_dist(correct.phyloseq)
  }else{
    correct.dist = distance(physeq = correct.phyloseq,method = method,...)
  }
  perm.dist = mclapply(1:dim(perm.GCN)[2],function(k){
    this.GCN = round(perm.GCN[,k])
    this.table = correct_abundance_table(abundance.table,this.GCN)
    this.phyloseq = phyloseq(otu_table(this.table,taxa_are_rows = FALSE))
    if(method=="Aitchison"){
      this.dist = calculate_Aitchison_dist(this.phyloseq)
    }else{
      this.dist = distance(physeq = this.phyloseq,method = method,...)
    }
    return(this.dist)
  },mc.cores = numCores)
  perm.PERMANOVA = mclapply(perm.dist,function(this.dist){
    adonis(formula = this.dist~V1,data = metadata,permutations = 999)
  },mc.cores = numCores)
  correct.PERMANOVA = adonis(formula = correct.dist~V1,data = metadata,permutations = 999)
  R2.CI = quantile(sapply(perm.PERMANOVA,function(this.PERMANOVA){
    this.PERMANOVA$aov.tab$R2[1]
  }),probs = c(alpha/2,1-alpha/2))
  CI.dist = sapply(1:length(perm.dist[[1]]),function(j){
    quantile(sapply(perm.dist,function(x){x[j]}),probs = c(alpha/2,1-alpha/2))
  })
  return(list(dist = correct.dist,CI=CI.dist,
              correct= correct.PERMANOVA,
              R2.CI = R2.CI))
}
#' @export
randomforest_cross_validation = function(predictor,response,strata=mod(1:length(response),10)){
  indice = lapply(unique(strata),function(x){which(strata==x)})
  # random forest with cross validation
  rfFit = lapply(indice,function(idx){
    res = numeric(length(response))
    mdl = randomForest::randomForest(x = predictor[-idx,], y = factor(response)[-idx], importance=TRUE, ntree = 500)
    res[idx] = predict(mdl, predictor[idx,])
    impt = mdl$importance[,'MeanDecreaseAccuracy']
    return(list(model=mdl,prediction=res,importance=impt))
  })
  # summarize prediction
  rfClasses = levels(factor(response))[apply(sapply(rfFit,function(x){x$prediction}),1,sum)]
  # find top predictors by feature importance
  rfImportance.mean = apply(sapply(rfFit,function(x){x$importance}),1,mean)
  rfImportance.sd = apply(sapply(rfFit,function(x){x$importance}),1,sd)
  top.predictor = sort(rfImportance.mean,decreasing = TRUE)
  return(list(confusion = table(rfClasses, response),rank=top.predictor,
              importance=list(mean=rfImportance.mean,sd=rfImportance.sd)))
}
#' @export
correct_randomforest_with_confidence = function(abundance.table,error,metadata,n=1000,numCores=1){
  perm.GCN = sapply(1:n,function(nn){
    this.GCN = sapply(error,function(x){rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)})
    this.GCN[this.GCN<=1] = 1
    return(this.GCN)
  })
  GCN = round(apply(perm.GCN,1,median))
  correct.table = correct_abundance_table(abundance.table,GCN)
  correct.rf = randomforest_cross_validation(predictor = as.matrix(correct.table),
                                             response = metadata$V1)
  perm.rf = mclapply(1:dim(perm.GCN)[2],function(k){
    this.GCN = round(perm.GCN[,k])
    this.table = abundance.table
    this.table = correct_abundance_table(abundance.table,this.GCN)
    randomforest_cross_validation(predictor = as.matrix(this.table),
                                  response = metadata$V1)
  },mc.cores = numCores)
  perm.rank = sapply(perm.rf,function(x){names(x$rank)})
  correct.support = sapply(1:length(correct.rf$rank),function(this.rank){
    local.support = sum(perm.rank[this.rank,]==names(correct.rf$rank)[this.rank])/n
    accumulative.support = sum(unlist(as.vector(perm.rank[1:this.rank,]))==names(correct.rf$rank)[this.rank])/n
    return(c(local=local.support,accumulate = accumulative.support))
  })
  return(list(correct= correct.rf,support = correct.support))
}
#' @export
correct_randomforest = function(abundance.table,GCN,metadata){
  correct_abundance_table(abundance.table,GCN)
  correct.rf = randomforest_cross_validation(predictor = as.matrix(correct.table),
                                             response = metadata$V1)
  return(correct= correct.rf)
}