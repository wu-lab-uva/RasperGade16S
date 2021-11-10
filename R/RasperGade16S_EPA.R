#' @title Run EPA from R terminal
#' @description Calls EPA-ng and insert query sequences into the reference
#' @export
#' @rdname insert_query_with_EPA
insert_query_with_EPA = function(seqs,tree,ref.seqs,model,save.path,numCores=0,max.out=99,intern=TRUE){
  if(missing(seqs)) stop("Path to trimmed query alignments must be provided!")
  if(missing(save.path)) save.path = "RasperGade16S_EPA/"
  if(!dir.exists(save.path)) dir.create(path = save.path,recursive = TRUE)
  if(missing(tree)) tree = system.file("extdata","reference.tre",package = "RasperGade16S",mustWork = TRUE)
  if(missing(ref.seqs)) ref.seqs = system.file("extdata","reference.trimmed.afa",package = "RasperGade16S",mustWork = TRUE)
  if(missing(model)) model = system.file("extdata","reference.RAxML.model.txt",package = "RasperGade16S",mustWork = TRUE)
  if(Sys.info()["sysname"]=="Windows") stop("EPA-ng currently not available on Windows.\nSee https://github.com/Pbdas/epa-ng for more information.")
  if(numCores<1){
    cmd = sprintf("epa-ng -t %s -s %s -q %s --model %s --preserve-rooting on --outdir %s --redo --filter-max %d",
                  tree,ref.seqs,seqs,model,save.path,max.out)
  }else{
    cmd = sprintf("epa-ng -t %s -s %s -q %s --model %s --preserve-rooting on --outdir %s --redo -T %d --filter-max %d",
                  tree,ref.seqs,seqs,model,save.path,numCores,max.out)
  }
  print(cmd)
  cmd.out=system(command = cmd,intern = intern)
  return(list(out=cmd.out,jplace=sprintf("%s/epa_result.jplace",save.path)))
}

#' @title Align and trim sequence
#' @description Align use HMMER 3
#' @export
#' @rdname align_with_HMM_and_trim
align_with_HMM_and_trim = function(seqs,hmm,mapali,mask,save.path){
  if(missing(seqs)) stop("Path to query sequences must be provided!")
  if(missing(save.path)) save.path = "RasperGade16S_align/"
  if(!dir.exists(save.path)) dir.create(path = save.path,recursive = TRUE)
  if(missing(hmm)) hmm = system.file("extdata","16S_core_gg.hmm",package = "RasperGade16S",mustWork = TRUE)
  if(missing(mapali)) mapali = system.file("extdata","aligned_remove_dot.fasta",package = "RasperGade16S",mustWork = TRUE)
  if(missing(mask)) mask = RasperGade16S.GG.13.8.mask.keys
  if(Sys.info()["sysname"]=="Windows") stop("HMMER3 currently not available on Windows.\nSee http://hmmer.org/ for more information.")
  cmd = sprintf("hmmalign --trim --dna -o %s/seq.align --mapali %s %s %s",
                save.path,mapali,hmm,seqs)
  print(cmd)
  cmd.out=system(command = cmd,intern = TRUE)
  cmd = sprintf("esl-reformat -o %s/seq.afa -u --gapsym=- afa %s/seq.align",save.path,save.path)
  print(cmd)
  cmd.out=system(command = cmd,intern = TRUE)
  if(is.null(mask)){
    file.copy(from = sprintf("%s/seq.afa",save.path),to = sprintf("%s/trimmed.afa",save.path),
              overwrite = TRUE)
  }else{
    trim_sequence_with_mask(align = sprintf("%s/seq.afa",save.path),
                            trimmed.align = sprintf("%s/trimmed.afa",save.path),
                            mask=mask)
  }
  return(list(out=cmd.out,afa = sprintf("%s/trimmed.afa",save.path)))
}

#' @title Trim 16S sequence based on GreenGene 13.8 16S mask
#' @description Trim 16S sequence based on GreenGene 13.8 16S mask
#' @export
#' @rdname trim_sequence_with_mask
trim_sequence_with_mask = function(align,trimmed.align,mask){
  if(missing(mask)){
    mask = RasperGade.GG.13.8.mask.keys
  }else{
    if(length(mask)<=1) mask = readRDS(mask)
    }
  afa = seqinr::read.fasta(file = align,seqtype = "DNA",as.string = FALSE,forceDNAtolower = FALSE)
  afa.id = unname(sapply(afa,function(x){attr(x,"name")}))
  matchedID = match(attr(mask,"seq_id"),afa.id)
  if(any(is.na(matchedID))) stop("Sequences in mapped alignment missing!\n")
  ref.afa = sapply(afa[matchedID],function(x){res=x;res[x=="."]="-";return(res)})
  pos.key = apply(ref.afa,1,function(x){toupper(paste0(x,collapse = ""))})
  inseq=seqinr::c2s(rep("-",length(matchedID)))
  pos2mask = which(!is.na(match(pos.key,mask)))
  afa = afa[-matchedID]
  trim.afa = lapply(afa,function(x){x[pos2mask]})
  seqinr::write.fasta(sequences = trim.afa,names = sapply(afa,function(x){attr(x,"name")}),
                      file.out = trimmed.align,nbchar = length(mask))
}
