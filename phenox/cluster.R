setwd('../data')
library(pvclust)
library(ape)
library(circlize)
library(ComplexHeatmap)

### read in the file
da <- read.csv('gds.gene.mat.csv',row.names = 1)
da <- da[,colSums(da) < 0.8*nrow(da)]

# Ward Hierarchical Clustering with Bootstrapped p values

fit <- pvclust(t(da), method.hclust="ward.D2",nboot=1000,
               method.dist="euclidean")
pdf('hcluster.pdf',height = 4.3, width = 5.8)
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)
dev.off()


### Output a tree file
hc <- fit$hclust
class(hc) # must be hclust class
my_tree <- as.phylo(hc) 
write.tree(phy=my_tree, file="tree.newick.txt")

### Generate a membership of the cluster
clu <- pvpick(fit)$cluster
names(clu) <- paste0('cluster',seq(1:length(clu)),"_")
uclu <- unlist(clu)
names(uclu) <- gsub('_[0-9]*$','',names(uclu))

memb <- data.frame(cluster=names(uclu),gdsid = uclu)
write.csv(memb,
          'memb.csv',row.names = F)

#### output heatmap 
type <- rep("Not in a Cluster",nrow(da))
type[match(memb$gdsid,rownames(da))] <- as.character(memb$cluster[match(memb$gdsid,rownames(da))])
type <- as.factor(type)

pdf('heatmap.pdf',height = 4.3, width = 5.8)
Heatmap(t(da), cluster_columns = hc, cluster_rows = F,
        show_row_names = FALSE,
        heatmap_legend_param = list(title = "Presence"))
dev.off()
