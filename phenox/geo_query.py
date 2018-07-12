#!/usr/bin/env python

'''
Functions for querying ncbi for geoprofiles data and gds-pubmed

See main, eventually tests will be added for this module
'''

import logging
import sys
import time
import Bio.Entrez as Entrez
try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2


# global variables
Entrez.email = A.N.Other@example.com #maybe ask for in argparse
Entrez.tool = phenotypeXpression
efetch_batch = 5000
elink_batch = 1000 # try it though, should be able to handle 10000
database = 'geoprofiles'

# timing compatibility with python2
def timing_tool():
    try:
        return time.perf_counter()
    except:
        return time.time()

# http attempts loop for biopython, assumes Bio.Entrez as Entrez
def http_attempts(query_type, **kwargs):
    for attempt in range(4): # HTTPError and retry three times
        try:
            fetch_handle = getattr(Entrez, query_type)(**kwargs)
        except ValueError as oops:
            logging.warning('Empty set: Likely total number = batch size')
            break
        except HTTPError as err:
            if 500 <= err.code <= 599:
                logging.info("Attempt {} of 4".format(attempt + 1))
                logging.warning(err)
                time.sleep(15)
            elif err.code == 400:
                logging.info("Attempt {} of 4".format(attempt + 1))
                logging.warning(err)
                sys.tracebacklimit = None
                time.sleep(15)
            else:
                raise
        else:
            break
    else:
        logging.critical('Failed to connect to NCBI 4 times exiting...')
        sys.exit()
    return fetch_handle

# epost list to e-utilities history server for target db.
# Can't elink over 1000, no max on epost
def post_ncbi(query_type, **kwargs):
    logging.debug('{} to ncbi using kwargs: {}'.format(query_type, kwargs))
    fetch_handle = http_attempts(query_type, **kwargs)
    query = Entrez.read(fetch_handle)
    fetch_handle.close()

    webenv = query['WebEnv']
    query_key = query['QueryKey']
    count = int(query['Count'])

    logging.debug('returned webenv: {} and query key: {}'
                  .format(webenv, query_key))
    return webenv, query_key, count
    
# Utilized for batch efetch/esummary from ncbi using history server.
def batch_ncbi(query_type, query_results, count1, **kwargs):
    # number in history list
    count = count1

    # sliding window of history list in batch size
    for start in range(0, count, efetch_batch):
        end = min(count, start + efetch_batch)
        logging.info("Going to download record {} to {}"
                     .format(start + 1, end))

        logging.debug('{} run with {}'
                      .format(query_type, dict(retstart=start, **kwargs)))
        fetch_handle = http_attempts(query_type,
                                     **dict(retstart=start, **kwargs))
        query_results.append(Entrez.read(fetch_handle))
        fetch_handle.close()
    return

# batch the queries locally, then epost->etuility.
# elink can't handle batch size > 1000, and doesn't support retstart,retmax.
# this is then the only way to epost -> elink,egquery,etc.
def batch_local(query_type, id_list, **kwargs):
    count = len(id_list) #non-redundant gds list
    result_list = []

    for start in range(0, count, elink_batch):
        end = min(count, start + elink_batch)
        logging.info("Looking up PMIDs for GDSs {} to {}"
                     .format(start + 1, end))
        
        # epost each batch from original db
        webenv1, query_key1 = post_ncbi('epost', db=kwargs['dbfrom'],
                id=",".join(id_list[start:end]))
        logging.debug('{} run with {}'.format(query_type, dict(
                webenv=webenv1, query_key=query_key1, **kwargs)))
        # elink each batch
        fetch_handle = http_attempts(query_type, **dict(
                webenv=webenv1, query_key=query_key1, **kwargs))
        record = Entrez.read(fetch_handle)
        fetch_handle.close()
        # for each Link'Id', dig down through dict add to list
        result_list.extend(
                [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']])
    return result_list

# getting xml variation files for query_results list, 
# START HERE
def get_ncbi_docsum(mesh_term, database, query_results):
    logging.debug('looking up id list -> {}'.format(id_list))
    start = timing_tool()
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
    return

