import geo_query
import logging
import sys
import time
import Bio.Entrez as Entrez
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


# global variables
Entrez.email = 'A.N.Other@example.com' #maybe ask for in argparse
Entrez.tool = phenotypeXpression
efetch_batch = 5000
elink_batch = 100 # try it though, should be able to handle 10000
database = 'geoprofiles'

def get_ncbi_docsum(mesh_term, database):
    logging.debug('looking up id list -> {}'.format(id_list))
    start = timing_tool()
    query_results = []
    staging_list = []

    # First, esearch using mesh_term, and keep results on ePost history server
    webenv1, query_key1, count1 = post_ncbi('esearch', db=database, usehistory='y',
            query='"up down genes"[filter] AND "{}"'.format(mesh_term))

    # Then, efetch in batches of efetch_batch
    batch_ncbi('efetch', staging_list, count1, db=database,
               rettype='docsum', retmax=efetch_batch, webenv=webenv1,
               query_key=query_key1)
    
    # append to passed list
    [query_results.append(batch) for batch in staging_list] 

    '''
    Now there needs to be some parsing to get the specific data out of the query_results
    list of dicts and it needs to grab PMC or PMIDs from a non-redundant set of gds's
    check on the 1:1 correlation, should be one for each. Deliver {gds: PMID} dict for
    term mining. Link page says gds_pubmed can support 10,000, but it may not (plus you
    may have more than 10k results). Use batch local above to split list of gds for
    elink lookup.
    
    '''
    # pmid_list = batch_local('elink', id_list, db='pubmed', dbfrom='gds', linkname='gds_pubmed')
    # verify that len(pmid_list) == len(gds_list)

    # stop timing and log
    stop = timing_tool()
    logging.info('Download time: {} min, Batches run -> {}'
                 .format(((stop-start)/60), len(query_results)))
    return query_results

def batch_local(query_type, id_list, **kwargs):
    count = len(id_list) #non-redundant gds list
    result_dict = {}

    for start in range(0, count, elink_batch):
        end = min(count, start + elink_batch)
        logging.info("Looking up PMIDs for GDSs {} to {}"
                     .format(start + 1, end))
        
        # elink each batch
        fetch_handle = http_attempts(query_type, **dict(
                id=','.join(id_list[start:end]), **kwargs))
        record = Entrez.read(fetch_handle)
        fetch_handle.close()
        # for each Link'Id', dig down through dict add to list
        for el in record:
            result_dict[el['IdList'][0]] = el['LinkSetDb'][0]['Link'][0]['Id']
    return result_dict

def gdsdict_from_profile(query_results):
    gdsdict = {}
    genefreq = {}
    for docsum in query_results['DocumentSummarySet']['DocumentSummary']:
        gdsid = docsum['GDS']
        gids = docsum['ENTREZ_GENE_ID'].split(";")
        #gid = docsum['ENTREZ_GENE_ID'].split(";")
        #print(gdsid, gids,docsum.attributes['uid'])
        for gid in gids:
            if len(gid) > 0:
                if gid in genefreq:
                    genefreq[gid] += 1
                else:
                    genefreq[gid] = 1
                if gdsid in gdsdict:
                    if gid in gdsdict[gdsid]:
                        gdsdict[gdsid][gid] += 1
                    else:
                        gdsdict[gdsid][gid] = 1
                else:
                    gdsdict[gdsid] ={}
                    if gid in gdsdict[gdsid]:
                        gdsdict[gdsid][gid] += 1
                    else:
                        gdsdict[gdsid][gid] = 1
    return gdsdict,genefreq

def genedict_from_profile(query_results):
    genedict = {}
    for docsum in query_results['DocumentSummarySet']['DocumentSummary']:
        gidlist = docsum['ENTREZ_GENE_ID'].split(";")
        gnamelist = docsum['geneName'].split("<:>")
        for (x,y) in zip(gidlist,gnamelist):
            genedict[x] = y
    return genedict

'''
Commands to convert gds ids to pubmed id, it will return a dictionary with the gds IDs as the keys, and pubmed id as the values.
'''
# batch_local('elink', gdslist, db='pubmed', dbfrom='gds', linkname='gds_pubmed')
    