from datetime import datetime

def meta_from_gds(self,gds_list: list) -> Dict:
    """
    Get meta information, such as gds submission time, n_samples, platform from gds.
    gds_list is from gdsdict.keys()
    - This can replace the get_pubmed_ids()
    concerns:
    >> do not know how to get rid of the IntergerElement
    >> now the publication time is in seconds, we can convert to days or years?
    """
    qout = Entrez.read(Entrez.esummary(db="gds", id=",".join(gds_list)))
    return {sm['Id']:[sm['n_samples'],sm['PubMedIds'],sm['GPL'],datetime.strptime(sm['PDAT'],'%Y/%m/%d').timestamp()] for sm in qout}
    

def gNdict_from_profile(self, query_results) -> Tuple:
    """
    Generate GDS data dict from profile data
    using <gene name> as dict key
    :param query_results:
    :return:
    """
    gds_dict = dict()
    gene_freq = defaultdict(int)
    for batch in query_results:
        for docsum in tqdm.tqdm(batch['DocumentSummarySet']['DocumentSummary'],desc="GEO Datasets"):
            gdsid = docsum['GDS']
            gnamelist = docsum['geneName'].split("<:>")
            for gid in gnamelist:
                if len(gid) > 0:
                    gene_freq[gid.upper()] += 1
                    if gdsid not in gds_dict:
                        gds_dict[gdsid] = defaultdict(int)
                    gds_dict[gdsid][gid.upper()] += 1
    return gds_dict, gene_freq
