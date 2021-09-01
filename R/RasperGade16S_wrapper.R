#' @title Predict 16S rRNA GCN from sequences
#' @description A wrapper function to "one-click-and-run" the default pipeline
#' @export
#' @rdname predict_16SGCN_from_sequences
predict_16SGCN_from_sequences = function(seqs){
  if(missing(seqs)){
    seqs = system.file("extdata/Demo","demo.SILVA.fasta",package="RasperGade16S",mustWork=TRUE)
    cat(sprintf("No FASTA file supplied, running prediction using demo sequences in\n%s\n\n",seqs))
    }
  align.out = align_with_HMM_and_trim(seqs=seqs)
  epa.out = insert_query_with_EPA(seqs="RasperGade16S_align/trimmed.afa")
  insert.locations = parse_jplace(epa.out$jplace,split = 1)
  insert.prediction = lapply(insert.locations,function(this.insert){
    this.res = predictHiddenStateWithPE(FMR = RasperGade16S.refdata$FMR,
                                        query.keys = this.insert$hash,laplace = FALSE)
    return(this.res)
  })
  insert.res = list(hsp=do.call(rbind,lapply(insert.prediction,function(x){x$hsp})),
                    error=do.call(c,lapply(insert.prediction,function(x){x$error})))
  insert.discrete.res = discretizeResult(res = insert.res$hsp,
                                         error = insert.res$error,laplace = FALSE)
  insert.GCN = insert.discrete.res$x
  names(insert.GCN) = insert.discrete.res$label
  return(list(tab=insert.discrete.res[,-1],GCN=insert.GCN,error = insert.res$error))
}
