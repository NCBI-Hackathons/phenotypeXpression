import sys
import logging
import time
import tqdm
from collections import defaultdict
import Bio.Entrez as Entrez
from urllib.error import HTTPError


# class for querying GEO databases
class GEOQuery:
    def __init__(self, email, tool="phenotypeXpression", efetch_batch=5000, elink_batch=100):
        Entrez.email = email
        Entrez.tool = tool
        self.efetch_batch = efetch_batch
        self.elink_batch = elink_batch
        self.db = 'geoprofiles'

    @staticmethod
    def timing_tool():
        return time.perf_counter()

    @staticmethod
    def http_attempts(query_type, **kwargs):
        """
        HTTP attempts
        :param query_type:
        :param kwargs:
        :return:
        """
        fetch_handle = None

        for attempt in range(4):  # HTTPError and retry three times
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

    def post_ncbi(self, query_type, **kwargs):
        """
        Post request to NCBI
        :param query_type:
        :param kwargs:
        :return:
        """
        logging.debug('{} to ncbi using kwargs: {}'.format(query_type, kwargs))
        fetch_handle = self.http_attempts(query_type, **kwargs)
        query = Entrez.read(fetch_handle)
        fetch_handle.close()

        webenv = query['WebEnv']
        query_key = query['QueryKey']
        count = int(query['Count'])

        logging.debug('returned webenv: {} and query key: {}'
                      .format(webenv, query_key))
        return webenv, query_key, count

    # Utilized for batch efetch/esummary from ncbi using history server.
    def batch_ncbi(self, query_type, query_results, count1, **kwargs):
        # number in history list
        count = count1

        # sliding window of history list in batch size
        for start in tqdm.tqdm(range(0, count, self.efetch_batch)):
            end = min(count, start + self.efetch_batch)
            logging.info("Going to download record {} to {}"
                         .format(start + 1, end))

            logging.debug('{} run with {}'
                          .format(query_type, dict(retstart=start, **kwargs)))
            fetch_handle = self.http_attempts(
                query_type, **dict(retstart=start, **kwargs)
            )
            query_results.append(Entrez.read(fetch_handle))
            fetch_handle.close()
        return

    def get_ncbi_docsum(self, mesh_term, database):
        """
        Get document summaries from NCBI database
        :param mesh_term:
        :param database:
        :return:
        """
        start = self.timing_tool()
        staging_list = []
        query_results = []

        # First, esearch using mesh_term, and keep results on ePost history server
        webenv1, query_key1, count1 = self.post_ncbi(
            'esearch', db=database, usehistory='y',
            term='"up down genes"[filter] AND "{}"'.format(mesh_term)
        )

        # Then, efetch in batches of efetch_batch
        self.batch_ncbi('efetch', staging_list, count1, db=database,
                        rettype='docsum', retmax=self.efetch_batch, webenv=webenv1,
                        query_key=query_key1)

        # append to passed list
        [query_results.append(batch) for batch in staging_list]

        # stop timing and log
        stop = self.timing_tool()
        logging.info('Download time: {} min, Batches run -> {}'
                     .format(((stop - start) / 60), len(query_results)))

        return query_results

    def batch_local(self, query_type, id_list, **kwargs):
        """
        Batch local
        :param query_type:
        :param id_list:
        :param kwargs:
        :return:
        """
        count = len(id_list)
        result_dict = {}

        for start in range(0, count, self.elink_batch):
            end = min(count, start + self.elink_batch)
            logging.info("Looking up PMIDs for GDSs {} to {}"
                         .format(start + 1, end))

            # elink each batch
            fetch_handle = self.http_attempts(
                query_type, **dict(id=','.join(id_list[start:end]), **kwargs)
            )
            record = Entrez.read(fetch_handle)
            fetch_handle.close()
            # for each Link'Id', dig down through dict add to list
            for el in record:
                result_dict[el['IdList'][0]] = el['LinkSetDb'][0]['Link'][0]['Id']
        return result_dict

    def gdsdict_from_profile(self, query_results):
        """
        Generate GDS data dict from profile data
        :param query_results:
        :return:
        """
        gds_dict = dict()
        gene_freq = defaultdict(int)

        for batch in query_results:
            for docsum in tqdm.tqdm(batch['DocumentSummarySet']['DocumentSummary']):
                gdsid = docsum['GDS']
                gids = docsum['ENTREZ_GENE_ID'].split(";")
                for gid in gids:
                    if len(gid) > 0:
                        gene_freq[gid] += 1
                        if gdsid not in gds_dict:
                            gds_dict[gdsid] = defaultdict(int)
                        gds_dict[gdsid][gid] += 1

        return gds_dict, gene_freq

    def genedict_from_profile(self, query_results):
        """
        Generate gene data dict from profile data
        :param query_results:
        :return:
        """
        gene_dict = dict()

        for batch in query_results:
            for docsum in batch['DocumentSummarySet']['DocumentSummary']:
                gid_list = docsum['ENTREZ_GENE_ID'].split(";")
                gname_list = docsum['geneName'].split("<:>")
                for (x, y) in zip(gid_list, gname_list):
                    gene_dict[x] = y

        return gene_dict