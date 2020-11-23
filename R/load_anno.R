#' Load annotation data from gene ontology database
#'
#' @param orgdb package with sepcific species from AnnotationDbi
#' @param keytype ID class
#' @param ont ontology "BP" or "MF" or "CC"
#' @param drop logical value, TRUE means drop IEA evidence code
#'
#' @return a list contains four lists
#' @export
#'
#' @examples
#' tmp <- load_anno(org.Hs.eg.db, "ENTREZID", "BP", TURE)
#'
load_anno <- function(orgdb, keytype, ont, drop) {

  #调用GOSemSim包里的godata函数加载注释数据
  aa <- GOSemSim::godata(orgdb, keytype, ont, computeIC = F)
  #取有关GO terms and genes or gene products(proteins)的数据
  geneanno <- aa@geneAnno
  #drop参数表示是否保留evidence "IEA"
  if (drop) {
    geneanno <- geneanno[geneanno$EVIDENCE != "IEA", ]
  }
  #取所有节点和proteins
  nodes <- unique(geneanno$GO)
  proteins <- unique(geneanno[, keytype])

  #每个GO term对应的蛋白数量，用于计算ICA
  gocount <- as.list(table(geneanno$GO))

  #每个蛋白的term list
  pro_annotations <- lapply(proteins, function(e) {
    unique(geneanno[geneanno[, keytype] == e, "GO"])
  })
  names(pro_annotations) <- proteins
  #以列表返回多个值
  res <- list(nodes, proteins, gocount, pro_annotations)
  return(res)
}

#加载注释数据，得到term和proteins的对应关系
#orgdb 表示所选物种，keytype表示所求类别
#ont表示ontology,drop关于"IEA"
