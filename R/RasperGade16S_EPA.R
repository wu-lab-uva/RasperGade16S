#' @title Run EPA from R terminal
#' @description 
#' @export
#' @rdname insert_query_with_EPA
insert_query_with_EPA = function(seqs,tree,ref.seqs,model,save.path,numCores=1,max.out=1){
  if(missing(save.path)) save.path = "RasperGade16S_EPA/"
  if(!dir.exists(save.path)) dir.create(path = )
  if(missing(tree))
  cmd = sprintf("epa-ng -t %s -s %s -q %s --model %s --preserve-rooting on --outdir %s --redo -T %d --filter-max %d",
                tree,ref.seqs,seqs,model,save.path,numCores,max.out)
  cmd.out=system(command = cmd,intern = TRUE)
  return(list(out=cmd.out,jplace=sprintf("%sepa_result.jplace",save.path)))
}

#' @title Align and trim query sequence
#' @description 
#' @export
#' @rdname align_and_trim_query
align_and_trim_query = function(seqs,tree,ref,model,save.paths,numCores=1,max.out=1){
  if(missing(save.path)) stop("Save path has to be provided!")
  cmd = sprintf("epa-ng -t %s -s %s -q %s --model %s --preserve-rooting on --outdir %s --redo -T %d --filter-max %d",
                tree,ref,seqs,model,save.path,numCores,max.out)
  cmd.out=system(command = cmd,intern = TRUE)
  return(list())
}