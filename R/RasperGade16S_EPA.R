#' @title Run EPA from R terminal
#' @description 
#' @export
#' @rdname insert_query_with_EPA
insert_query_with_EPA = function(seqs,tree,ref.seqs,model,save.path,numCores=1,max.out=1,intern=TRUE){
  if(missing(seqs)) stop("Path to query sequences must be provided!")
  if(missing(save.path)) save.path = "RasperGade16S_EPA/"
  if(!dir.exists(save.path)) dir.create(path = save.path)
  if(missing(tree)) tree = system.file("extdata","reference.tre",package = "RasperGade16S",mustWork = TRUE)
  if(missing(ref.seqs)) ref.seqs = system.file("extdata","reference.trimmed.afa",package = "RasperGade16S",mustWork = TRUE)
  if(missing(model)) model = system.file("extdata","reference.RAxML.model.txt",package = "RasperGade16S",mustWork = TRUE)
  if(Sys.info()["sysname"]=="Windows") stop("EPA-ng currently not available on Windows.\nSee https://github.com/Pbdas/epa-ng for more information.")
  cmd = sprintf("epa-ng -t %s -s %s -q %s --model %s --preserve-rooting on --outdir %s --redo -T %d --filter-max %d",
                tree,ref.seqs,seqs,model,save.path,numCores,max.out)
  cmd.out=system(command = cmd,intern = intern)
  return(list(out=cmd.out,jplace=sprintf("%s/epa_result.jplace",save.path)))
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