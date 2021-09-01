#' Default reference tree for RasperGade16S
#'
#' A phylo-class object containing 6408 reference tips/genomes
#'
#' @format A phylo-class object with 6408 tips, binary and rooted:
#' \describe{
#'   \item{tip.label}{reference genome names, in NCBI RefSeq accession number}
#'   \item{edge.length}{phylogenetic distance between nodes, in substitutions/site}
#'   \item{edge}{parent node and descendant node of an edge}
#'   \item{Nnode}{the number of nodes}
#' }
#' @source \url{https://github.com/wu-lab-uva/16S-rRNA-GCN-Predcition}
"RasperGade16S.reftree"

#' Default reference GCN for RasperGade16S
#'
#' A named vector of 6408 elements
#'
#' @format A named numeric vector with 6408 elements:
#' \describe{
#'   \item{values}{16S rRNA GCN of reference genomes}
#'   \item{names}{reference genome names, in NCBI RefSeq accession number}
#' }
#' @source \url{https://github.com/wu-lab-uva/16S-rRNA-GCN-Predcition}
"RasperGade16S.refGCN"

#' Predicted SILVA 132 GCN (Bacteria only) by RasperGade16S
#'
#' A named vector of 592605 elements
#'
#' @format A named numeric vector with 592605 elements:
#' \describe{
#'   \item{values}{16S rRNA GCN of SILVA bacterial OTUs}
#'   \item{names}{reference genome names, in NCBI RefSeq accession number}
#' }
#' @source \url{https://github.com/wu-lab-uva/16S-rRNA-GCN-Predcition}
"RasperGade16S.SILVA.GCN"

