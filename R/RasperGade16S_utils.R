#' @title Convert absolute abundance to relative abundance
#' @description 
#' @export
#' @rdname relative_abundance
relative_abundance = function(abundance){
  return(abundance/sum(abundance))
}

#' @export
sample_GCN_from_error_distribution = function(error,lower.edge="bounce"){
  this.GCN = round(sapply(error,function(x){
    rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)
  }))
  if(lower.edge=="bounce"){
    this.GCN[this.GCN<1] = 2-this.GCN[this.GCN<1]
  }else{
    this.GCN[this.GCN<1] = 1
  }
  return(this.GCN)
}

#' @export
correct_GCN = function(gene.abundance,GCN,normalize =TRUE){
  cell.abundance = gene.abundance/GCN
  if(normalize) cell.abundance = relative_abundance(cell.abundance)*sum(gene.abundance)
  return(cell.abundance)
}

#' @export
correct_abundance_table = function(abundance,GCN,normalize=TRUE){
  as.data.frame(t(apply(abundance,1,correct_GCN,GCN = GCN,normalize=normalize)))
}

#' @export
correct_abundance_table_with_SILVA = function(abundance,normalize=TRUE){
  as.data.frame(t(apply(abundance,1,correct_GCN,GCN = GCN,normalize=normalize)))
}