#' @title Run EPA from R terminal
#' @description Calls EPA-ng and insert query sequences into the reference
#' @export
#' @rdname insert_query_with_EPA
insert_query_with_EPA = function(seqs,tree,ref.seqs,model,save.path,numCores=1,max.out=1,intern=TRUE){
  if(missing(seqs)) stop("Path to trimmed query alignments must be provided!")
  if(missing(save.path)) save.path = "RasperGade16S_EPA/"
  if(!dir.exists(save.path)) dir.create(path = save.path,recursive = TRUE)
  if(missing(tree)) tree = system.file("extdata","reference.tre",package = "RasperGade16S",mustWork = TRUE)
  if(missing(ref.seqs)) ref.seqs = system.file("extdata","reference.trimmed.afa",package = "RasperGade16S",mustWork = TRUE)
  if(missing(model)) model = system.file("extdata","reference.RAxML.model.txt",package = "RasperGade16S",mustWork = TRUE)
  if(Sys.info()["sysname"]=="Windows") stop("EPA-ng currently not available on Windows.\nSee https://github.com/Pbdas/epa-ng for more information.")
  cmd = sprintf("epa-ng -t %s -s %s -q %s --model %s --preserve-rooting on --outdir %s --redo -T %d --filter-max %d",
                tree,ref.seqs,seqs,model,save.path,numCores,max.out)
  cmd.out=system(command = cmd,intern = intern)
  return(list(out=cmd.out,jplace=sprintf("%s/epa_result.jplace",save.path)))
}

#' @title Align and trim sequence
#' @description Align use HMMER 3
#' @export
#' @rdname align_with_HMM_and_trim
align_with_HMM_and_trim = function(seqs,hmm,mapali,save.path){
  if(missing(seqs)) stop("Path to query sequences must be provided!")
  if(missing(save.path)) save.path = "RasperGade16S_align/"
  if(!dir.exists(save.path)) dir.create(path = save.path,recursive = TRUE)
  if(missing(hmm)) hmm = system.file("extdata","16S_core_gg.hmm",package = "RasperGade16S",mustWork = TRUE)
  if(missing(mapali)) mapali = system.file("extdata","aligned_remove_dot.fasta",package = "RasperGade16S",mustWork = TRUE)
  if(Sys.info()["sysname"]=="Windows") stop("HMMER3 currently not available on Windows.\nSee https://http://hmmer.org/ for more information.")
  cmd = sprintf("hmmalign --trim --dna -o %s/seq.align --mapali %s %s %s",
                save.path,mapali,hmm,seqs)
  cmd.out=system(command = cmd,intern = TRUE)
  cmd = sprint("esl-reformat -o %s/seq.afa -u --gapsym=- afa %s/seq.align",save.path,save.path)
  trim_sequence_with_mask(align = sprintf("%s/seq.afa",save.path),trimmed.align = sprintf("%s/trimmed.afa",save.path))
  return(list(out=cmd.out,afa = sprintf("%s/trimmed.afa",save.path)))
}