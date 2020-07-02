packages = c(
  "BiocManager",
  'ggbeeswarm'
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
}

if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways")
}
library(rWikiPathways)
pathwaynames.hs <- listPathwayNames(organism = "Homo sapiens")
pathwaynames.hs[grepl(pattern = "MET", x = pathwaynames.hs)]
