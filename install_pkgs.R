packages = c(
  "BiocManager",
  'ggbeeswarm'
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
}
