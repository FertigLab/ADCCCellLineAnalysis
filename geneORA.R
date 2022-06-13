require(goSTAG)
require(GO.db)

myPGE = function (gene_lists, go_terms, filter_method = "pval", 
          significance_threshold = 0.05, p.adjust_method = "BH") 
{
  if (!is.null(go_terms[["ALL"]])) {
    all_genes = go_terms[["ALL"]]
    go_terms[["ALL"]] = NULL
  }
  else {
    stop("All annotated genes entry missing in GO terms list")
  }
  gene_list_sizes = vapply(gene_lists, function(x) {
    length(intersect(toupper(x), toupper(all_genes)))
  }, integer(1))
  universe_size = length(unique(toupper(all_genes)))
  pvals = matrix(nrow = length(go_terms), ncol = length(gene_lists))
  for (i in seq_len(length(go_terms))) {
    m = length(intersect(toupper(go_terms[[i]]), toupper(all_genes)))
    pvals[i, ] = vapply(seq_len(length(gene_lists)), function(j) {
      x = length(intersect(toupper(go_terms[[i]]), toupper(gene_lists[[j]])))
      k = gene_list_sizes[j]
      N = universe_size
      phyper(x - 1, m, N - m, k, lower.tail = FALSE)
    }, numeric(1))
  }
  
  adjusted_pvals = matrix(nrow = length(go_terms), ncol = length(gene_lists))
  for (i in seq_len(length(gene_lists))) {
    adjusted_pvals[, i] = p.adjust(pvals[, i], method = p.adjust_method)
  }
  if (toupper(filter_method) == "PVAL") {
    final_matrix = pvals[apply(pvals < significance_threshold, 
                               1, sum) > 0, ]
    final_matrix[final_matrix >= significance_threshold] = 1
  }
  else if (toupper(filter_method) == "P.ADJUST") {
    final_matrix = as.vector(adjusted_pvals)
    names(final_matrix) = names(go_terms)
    final_matrix[adjusted_pvals >= significance_threshold] = 1
    final_matrix = final_matrix[apply(adjusted_pvals < significance_threshold, 
                                      1, sum) > 0]
  }
  else {
    stop("Invalid filter method")
  }
  
  return(final_matrix)
}

geneORA = function(geneList, significance_threshold = 0.01, species = "human", p.adjust_method = "BH") {
  stagList = list(geneList = geneList)
  go_terms = loadGOTerms(species = species)
  enrichment_matrix = myPGE(stagList, go_terms, filter_method = "p.adjust",
                                          significance_threshold = significance_threshold,
                                          p.adjust_method = p.adjust_method)
  GOTerms <- as.list(GOTERM)
  goIDs = names(enrichment_matrix)
  terms = unlist(lapply(goIDs, function(x){Term(x)}))
  names(enrichment_matrix) = terms
  return(enrichment_matrix)
}
  

