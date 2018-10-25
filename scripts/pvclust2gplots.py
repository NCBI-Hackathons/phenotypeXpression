def call_r_clustering(self, gds_py):
    """
    Call R script
    :param gds: dataframe with GDS data
    :return:
    """
    # load all R libraries
    base = importr("base")
    pvclust = importr("pvclust")
    graphics = importr("graphics")
    grdevices = importr("grDevices")
    ape = importr("ape")
    gplots = importr("gplots")

    # convert pandas df to R df
    with localconverter(default_converter + pandas2ri.converter) as cv:
        gds = pandas2ri.py2ri(gds_py)

    matrix = base.as_matrix(gds)
    mat_trans = matrix.transpose()

    # cluster over studies
    print("Clustering on Studies...")
    fit = pvclust.pvclust(mat_trans, nboot=1000, method_hclust="ward.D2", method_dist="euclidean")

    # write clustering output to pdf
    grdevices.pdf(self.hcluster_file, paper="a4")
    graphics.plot(fit)
    pvclust.pvrect(fit, alpha=.95)
    grdevices.dev_off()
    print("Clustering diagram written to {}".format(self.hcluster_file))

    # write heatmap with pvclust to pdf 
    grdevices.pdf(self.heatmap_file, paper="a4")
    gplots.heatmap_2(mat_trans,Rowv = fit,dendrogram = "col",col=gplots.bluered(100), labRow = False, trace="none")
    grdevices.dev_off()
    print("Heatmap written to {}".format(self.heatmap_file))

    
    # write cluster tree output to tree file
    hc = fit.rx2("hclust")
    tree = ape.as_phylo(hc)
    ape.write_tree(phy=tree, file=self.tree_file)
    print("Cluster tree written to {}".format(self.tree_file))

    # extract cluster membership
    pvp = pvclust.pvpick(fit)
    clusters = pvp.rx2("clusters")

    cluster_members = defaultdict(list)

    if clusters[0] == rinterface.NULL:
        print("WARNING: Only one cluster!")
        cluster_members['cluster0'] = list(gds_py.index)
    else:
        for i, clust in enumerate(clusters):
            clust_name = "cluster{}".format(i)
            for pmid in clust:
                cluster_members[clust_name].append(pmid)

    # generate heatmap
    # self._generate_heatmap(gds_py, cluster_members)

    # generate distance graph
    self._generate_dist_graph(gds_py)

    return cluster_members