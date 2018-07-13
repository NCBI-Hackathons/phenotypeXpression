#### FUNCTION ####
installLibrary <- function(pkg_name){
  #### Check required package and install if the package is not in the library 
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg_name, repos='http://cran.us.r-project.org', quiet=TRUE)
  }
  require(pkg_name, character.only=TRUE, quietly=TRUE)
}

installBioLib <- function(pkg_name){
  source("https://bioconductor.org/biocLite.R")
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    biocLite(pkg_name)
  }
  require(pkg_name, character.only=TRUE, quietly=TRUE)
}
#### Required packages ######
.cran_packages <- c("pvclust", "ape", "optparse", "igraph")
.bio_packages <- c("circlize", "ComplexHeatmap")
sapply(c(.cran_packages), installLibrary)
sapply(c(.bio_packages), installBioLib)

#### Command line options #####
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="the directory where the input data file is.", metavar="character")
  )

opt_parser = OptionParser(usage = "usage: Rscript %prog -d <filedir>",
                          option_list=option_list,
                          description = "\t\t##### <<<< Graph generating R script for phenoXpress >>>> #####")
opt = parse_args(opt_parser)

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Please supply the path of the data directory as -d/--dir <filedir>.\n", call.=FALSE)
}
##### FUNCTION END #####

### read in the file  ####
da <- read.csv(paste0(opt$dir,'/gds.gene.mat.csv'),row.names = 1)
da <- da[,colSums(da) < 0.8*nrow(da)]

##### Ward Hierarchical Clustering with Bootstrapped p values  ####

fit <- pvclust(t(da), method.hclust="ward.D2",nboot=1000,
               method.dist="euclidean")
pdf(paste0(opt$dir,'/hcluster.pdf'),height = 4.3, width = 5.8)
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)
dev.off()

### Output a tree file  ####
hc <- fit$hclust
class(hc) # must be hclust class
my_tree <- as.phylo(hc) 
write.tree(phy=my_tree, file=paste0(opt$dir,"/tree.newick.txt"))

### Generate a membership of the cluster  ####
clu <- pvpick(fit)$cluster
names(clu) <- paste0('cluster',seq(1:length(clu)),"_")
uclu <- unlist(clu)
names(uclu) <- gsub('_[0-9]*$','',names(uclu))

memb <- data.frame(cluster=names(uclu),gdsid = uclu)
write.csv(memb,
          paste0(opt$dir,'/memb.csv'),row.names = F)

#### output heatmap ####
type <- rep("Not in a Cluster",nrow(da))
type[match(memb$gdsid,rownames(da))] <- as.character(memb$cluster[match(memb$gdsid,rownames(da))])
type <- as.factor(type)

pdf(paste0(opt$dir,'/heatmap.pdf'),height = 4.3, width = 5.8)
Heatmap(t(da), cluster_columns = hc, cluster_rows = F,
        show_row_names = FALSE,
        heatmap_legend_param = list(title = "Presence"))
dev.off()

#### gds distance graph ####
gds_dist = dist(da,method="binary")
pdf(file=paste0(opt$dir,'/GDS_dist_graph.pdf'), useDingbats=FALSE,pointsize=10)
g = graph.full(nrow(da))
V(g)$label = paste("GDS",rownames(da))
layout = layout.mds(g, dist = as.matrix(gds_dist))
plot(g, layout = layout, vertex.size = 1)
dev.off()