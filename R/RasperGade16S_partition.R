#' @export
#' @rdname get_connected_nodes
get_connected_nodes = function(tree,node,group){
  if(node<=Ntip(tree)) stop("Tips are not expected!")
  this.cluster = node
  this.outline = setdiff(getConnected(tree,node),1:Ntip(tree))
  while(length(this.outline)>0){
    this.outline = this.outline[(group[this.outline]==group[node])&
                                  (!this.outline%in%this.cluster)]
    this.cluster = c(this.cluster,this.outline)
    this.outline = do.call(c,lapply(this.outline,function(x){
      setdiff(getConnected(phy=tree,x = x),1:Ntip(tree))}))
  }
  num.transition = sum(sapply(this.cluster,function(x){
    this.parent = getLatestAncestor(phy = tree,x = x)
    if(is.na(this.parent)) this.parent = x
    this.des = getNextDescendants(phy = tree,x = x)
    this.trans = as.numeric(!this.parent%in%this.cluster) +
      sum(as.numeric(!this.des%in%c(this.cluster,1:Ntip(tree))))
  }))
  return(list(cluster=this.cluster,num.transition = num.transition,
              group = group[node],size=length(this.cluster)))
}

#' @export
#' @rdname check_partition
check_partition = function(tree,dAIC,group){
  remaining.nodes = unique(reorder.phylo(x = tree,order = "postorder",
                                         index.only = FALSE)$edge[,1])
  clusters = list()
  numCluster = 1
  while(length(remaining.nodes>0)){
    clusters[[numCluster]] = get_connected_nodes(tree = tree,node = remaining.nodes[1],
                                                group = group)
    clusters[[numCluster]]$dAIC = sum(dAIC[clusters[[numCluster]]$cluster-Ntip(tree)])*
      group[remaining.nodes[1]]+
      2*clusters[[numCluster]]$num.transition
    remaining.nodes = setdiff(remaining.nodes,clusters[[numCluster]]$cluster)
    numCluster = numCluster + 1
  }
  return(clusters)
}

#' @export
#' @rdname partition_phylogeny
partition_phylogeny_by_AIC_single = function(tree,dAIC,flip=FALSE){
  postorder.node = unique(reorder.phylo(x = tree,order = "postorder",
                                        index.only = FALSE)$edge[,1])
  all.group = numeric(Ntip(tree)+Nnode(tree))+(2*(as.numeric(dAIC[1]<0))-1)
  baseline.group = all.group
  #
  for(node in rev(postorder.node)){
    node.only.idx = node - Ntip(tree)
    this.des = getNextDescendants(tree,node)
    null.model = baseline.group[node]
    if((dAIC[node.only.idx]*null.model)>2){
      all.group[node] = -null.model
      all.group[this.des] = -null.model
      baseline.group[this.des] = -null.model
    }else{
      all.group[node] = null.model
      all.group[this.des] = null.model
      baseline.group[this.des] = null.model
    }
  }
  #
  if(flip){
    des.nodes = get_descendant_nodes_for_each_node(tree)
    for(ii in 1:length(postorder.node)){
      current.AIC = -sum(dAIC[all.group[Ntip(tree)+(1:Nnode(tree))]<0])+
        2*sum(all.group!=baseline.group)
      this.node = rev(postorder.node)[ii]
      this.node.only.idx = this.node - Ntip(tree)
      this.this.des = getNextDescendants(tree,this.node)
      new.group = all.group
      new.baseline.group = baseline.group
      new.group[this.node] = -all.group[this.node]
      new.baseline.group[this.this.des] = -all.group[this.node]
      for(node in rev(postorder.node[postorder.node%in%des.nodes[[this.node]]])[-1]){
        node.only.idx = node - Ntip(tree)
        this.des = getNextDescendants(tree,node)
        null.model = new.baseline.group[node]
        if((dAIC[node.only.idx]*null.model)>2){
          new.group[node] = -null.model
          new.group[this.des] = -null.model
          new.baseline.group[this.des] = -null.model
        }else{
          new.group[node] = null.model
          new.group[this.des] = null.model
          new.baseline.group[this.des] = null.model
        }
      }
      new.AIC = -sum(dAIC[new.group[Ntip(tree)+(1:Nnode(tree))]<0])+
        2*sum(new.group!=new.baseline.group)
      if(new.AIC<current.AIC){
        all.group = new.group
        baseline.group = new.baseline.group
      }
    }
  }
  #
  num.transition = sum(all.group!=baseline.group)
  return(list(group=all.group,
              num.transition = num.transition,
              dAIC=-sum(dAIC[all.group[Ntip(tree)+(1:Nnode(tree))]<0])+2*num.transition))
}

#' @export
#' @rdname get_partition_AIC
get_partition_AIC = function(tree,trait,models,group){
  reg.ancs = 
    reconstructAncestralStates(phy = tree,x = trait[tree$tip.label],
                               rate = total.process.variance(models[[1]]),
                               epsilon = models[[1]]["epsilon"])
  slow.ancs = 
    reconstructAncestralStates(phy = tree,x = trait[tree$tip.label],
                               rate = total.process.variance(models[[2]]),
                               epsilon = models[[2]]["epsilon"])
  reg.AIC = sapply(1:Nnode(tree),function(i){
    -2*log(dPEpoisnorm(x = reg.ancs$contrast[i],t = reg.ancs$l[i],
                       lambda=unname(models[[1]][1]),
                       size=unname(models[[1]][2]),
                       sigma=unname(models[[1]][3]),
                       epsilon=unname(reg.ancs$epsilon[i])))
  })
  slow.AIC = sapply(1:Nnode(tree),function(i){
    -2*log(dPEpoisnorm(x = slow.ancs$contrast[i],t = slow.ancs$l[i],
                       lambda=unname(models[[2]][1]),
                       size=unname(models[[2]][2]),
                       sigma=unname(models[[2]][3]),
                       epsilon=unname(slow.ancs$epsilon[i])))
  })
  postorder.node = unique(reorder.phylo(x = tree,order = "postorder",
                                        index.only = FALSE)$edge[,1])
  num.transition = sum(sapply(rev(postorder.node)[-1],function(x){
    group[getLatestAncestor(tree,x)]!=group[x]
  }))
  partition.AIC = sum(reg.AIC[group[Ntip(tree)+(1:Nnode(tree))]>0])+
    sum(slow.AIC[group[Ntip(tree)+(1:Nnode(tree))]<0])+
    2*num.transition
  return(partition.AIC)
}

#' @export
#' @rdname get_partition_transition
get_partition_transition = function(tree,group){
  postorder.node = unique(reorder.phylo(x = tree,order = "postorder",
                                        index.only = FALSE)$edge[,1])
  num.transition = sum(sapply(rev(postorder.node)[-1],function(x){
    group[getLatestAncestor(tree,x)]!=group[x]
  }))
  return(num.transition)
}

#' @export
#' @rdname refine_partition
refine_partition = function(tree,dAIC,group){
  current.partition = check_partition(tree,dAIC,group)
  current.dAIC = sapply(current.partition,function(x){x$dAIC})
  current.group = sapply(current.partition,function(x){x$group})
  last.pos.count = sum(current.dAIC>0)
  by.group=TRUE
  while(any(current.dAIC>0)){
    cat(sprintf("%d positive clusters remaining\n",last.pos.count))
    this.cluster = current.partition[[which.max(current.dAIC)]]
    this.group = this.cluster$group
    these.clusters = intersect(which(current.dAIC>0),which(current.group==this.group))
    if(!by.group) these.clusters = which.max(current.dAIC)
    group[do.call(c,lapply(current.partition[these.clusters],
                           function(x){x$cluster}))] = -this.group
    group[1:Ntip(tree)] = sapply(1:Ntip(tree),function(i){
      group[getLatestAncestor(phy = tree,x = i)]
    })
    current.partition = check_partition(tree,dAIC,group)
    current.dAIC = sapply(current.partition,function(x){x$dAIC})
    if(sum(current.dAIC>0)==last.pos.count) by.group = FALSE
    last.pos.count = sum(current.dAIC>0)
  }
  return(list(group=group,partition=current.partition))
}
#' @export
#' @rdname rescale_tree_by_partition_and_model
rescale_tree_by_partition_and_model = 
  function(tree,group,models,epsilon.as.PE = FALSE,min.value=1e-6){
    base.rate = total.process.variance(models[[1]])
    scale.coef = sapply(models,function(x){
      base.rate/total.process.variance(x)
    })
    for(i in 1:length(scale.coef)){
      tree$edge.length[group[tree$edge[,1]]==i] = 
        tree$edge.length[group[tree$edge[,1]]==i]/scale.coef[i]
    }
    new.epsilon = sapply(1:Ntip(tree),function(j){
      models[[group[j]]]["epsilon"]/2
    })
    new.model = models[[1]]
    new.model["epsilon"] = min.value
    scale.branch = scale.coef[group]
    epsilon.branch = rep(0,Ntip(tree)+Nnode(tree))
    if(epsilon.as.PE){
      for(i in 1:Ntip(tree)){
        tree$edge.length[tree$edge[,2]==i] =
          tree$edge.length[tree$edge[,2]==i]+new.epsilon[i]/base.rate
      }
      new.epsilon = rep(min.value,Ntip(tree))
      epsilon.branch = sapply(group,function(x){
        models[[x]]["epsilon"]/2/base.rate
      })
    }
    return(list(phy=tree,model= new.model,epsilon=unname(new.epsilon),
                scale.branch = scale.branch,epsilon.branch=epsilon.branch))
  }
#