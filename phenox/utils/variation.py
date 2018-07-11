#!/usr/bin/env python

'''
Below are the functions for the Variation file object which contains
variant info for metrics.

Clinotator - Clinical interpretation of ambiguous ClinVar annotations
Copyright (C) 2017  Robert R Butler III

See main, eventually tests will be added for this module
'''

import pprint
import logging
import datetime
import warnings
import decimal
from math import fsum
import numpy as np
import xml.etree.ElementTree as ET
import global_vars as g


# work through batch queries from ncbi and build list of variant objects
def query_parsing(variant_objects, query_results):
    
    for batch_index, batch in enumerate(query_results):
        clinvarresult = ET.fromstring(batch)
        for var_index, variationreport in enumerate(clinvarresult):
            variant = VariationClass(variationreport)
            variant_objects.append(variant)
        logging.debug('{} variant objects parsed in batch {}'
                      .format(var_index + 1, batch_index + 1))

# calculate the age of each individual assertion
def calculate_age(datestring):
    born = datetime.datetime.strptime(datestring, '%Y-%m-%d')
    today = datetime.date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month,
                                                                 born.day))

# calculate age weight of each assertion
def age_weight(age):

    if age < 2:
        return 1
    elif 2 <= age < 7:
        return ((11 - age) / 10)
    elif age >= 7:
        return 0.5

# when using a list of ages, catch empty list and return average age
def average_list_age(vid, age_list):
    #numpy RuntimeWarning for empty age_list, suppress warning and log 
    with warnings.catch_warnings():
        warnings.simplefilter('error', category=RuntimeWarning)

        try:
            return np.nanmean(age_list)

        except:
            logging.warning('VID: {} does not have valid clinical assertions!'
                          .format(vid))
            return None

# try a key to see if it is in the dict or warn
def key_test(test_dict, test_key):
    
    try:
        return test_dict[test_key]
    except KeyError:
        logging.warning('{} is not an expected key'.format(test_key))

# evaluates reclassification recommendation for CTRR
def reclassification_tree(ctps_index, cvcs_index):    
    index_diff = abs(ctps_index - cvcs_index)
    
    if index_diff == 0:
        return 0
    elif index_diff == 1 and cvcs_index in (0,1,2,4,5,6) and ctps_index != 3:
        return 1
    elif index_diff == 1 and cvcs_index in (0,1,2,4,5,6) and ctps_index == 3:
        return 2
    elif index_diff == 1 and cvcs_index == 3:
        return 2
    elif (index_diff == 2 and cvcs_index in (0, 1, 2)
            and ctps_index in (0, 1, 2)):
        return 2
    elif (index_diff == 2 and cvcs_index in (0, 1, 2)
            and ctps_index not in (0, 1, 2)):
        return 3
    elif index_diff == 2 and cvcs_index == 3:
        return 3
    elif (index_diff == 2 and cvcs_index in (4, 5, 6)
            and ctps_index in (4, 5, 6)):
        return 2
    elif (index_diff == 2 and cvcs_index in (4, 5, 6)
            and ctps_index not in (4, 5, 6)):
        return 3
    elif index_diff > 2:
        return 3
                

class VariationClass:
    
    def __init__(self, variationreport):
        self.VID = variationreport.attrib['VariationID']
        self.CVVT = variationreport.attrib['VariationType']
        self.allele_parse(variationreport)
        self.observation_parse(variationreport)
        self.assertion_table_stats(variationreport)
        self.analysis_stats()

    # parse the Allele subtree of variation report
    def allele_parse(self, variationreport):
        RS = []
        Alt = []
        vcf_match = []
        
        for index, alleles in enumerate(variationreport.findall('./Allele')):
            
            try:
                RS.append(alleles.find('./XRefList/XRef[@DB="dbSNP"]')
                          .get('ID'))
            except:
                RS.append('.')
            
            try:
                Alt.append(alleles
                           .find('./SequenceLocation[@Assembly="GRCh38"]')
                           .get('alternateAllele'))
            except:
                Alt.append('.')
                
            vcf_match.append('{}|{}'.format(RS[index], Alt[index]))
            
        logging.debug('RS list -> {}\nAlt list -> {}'.format(RS, Alt))
        self.rsID = RS
        self.CVAL = Alt
        
        if len(vcf_match) > 1:
            self.vcfmatch = vcf_match
        else:
            self.vcfmatch = '.'

    # parse the ObservationList subtree of variation report
    def observation_parse(self, variationreport):
        run_already = False

        for observation in variationreport \
                .findall('./ObservationList/Observation'):
            
            if (observation.get('VariationID') == self.VID and
                    not run_already):
                run_already = True
                reviewstat = observation.find('ReviewStatus').text
                self.CVSZ = g.star_dict[reviewstat]
                self.CVCS = observation \
                    .find('./ClinicalSignificance/Description').text
                self.CVLE = observation.find('./ClinicalSignificance') \
                    .get('DateLastEvaluated')

            elif (observation.get('VariationID') == self.VID and
                    run_already):
                logging.warning('{} has multiple observation fields in its re'
                             'cord omitting as an annotation error. Check rsi'
                             'd(s) {} manually'.format(self.VID, self.rsID))
                continue

            else:
                continue
        
    # parse the PhenotypeList subtree of variation report
    def pheno_parse(self, assertion, sig_key):
        pheno_list = []
        
        for phenotype in assertion.findall('./PhenotypeList/Phenotype'):
            pheno_list.append('{}({})'.format(phenotype.get('Name'), sig_key))
        
        return pheno_list

    # parse the ClinicalAssertionList subtree of variation report
    def assertion_table_stats(self, variationreport):
        raw_score = []
        age_list = []
        cvds_list = []

        for assertion in variationreport \
                .findall('./ClinicalAssertionList/GermlineList/Germline'):
            revstat_key = assertion.find('ReviewStatus').text
            score = key_test(g.cutoff, revstat_key)
            sigval_key = assertion.find('./ClinicalSignificance/Description') \
                .text
            sig_value = key_test(g.significance, sigval_key)

            if score > 0 and sig_value[0] != 0:
                try:
                    age = calculate_age(assertion
                                        .find('./ClinicalSignificance')
                                        .get('DateLastEvaluated'))
                except:
                    logging.warning('{} has a missing assertion date!'
                                 .format(self.VID))
                    continue
                    
                age_list.append(age)
                D = decimal.Decimal
                raw_score.append(float(D(str(score)) * D(str(sig_value[0]))
                                 * D(str(age_weight(age)))))
                logging.debug('score: {} sig_value: {} age_weight: {} age: {}'
                              .format(score, sig_value[0], age_weight(age), age))

                cvds_list += self.pheno_parse(assertion, sig_value[1])
            
        self.CVDS = ';'.join(cvds_list)
        if not self.CVDS:
            self.CVDS = '.'
        
        self.CVNA = len(age_list)
        self.CTAA = average_list_age(self.VID, age_list)
        # logging.debug('VID: {} -> age list size: {}, raw_score size: {}'
        #               .format(self.VID, len(age_list), len(raw_score)))

        if len(raw_score) >= 2:
            # logging.debug('raw_score list {}'.format(raw_score))
            self.CTRS = fsum(raw_score)
        else:
            self.CTRS = None
            logging.debug('VID: {} has fewer than 2 complete assertions, no C'
                          'TRS will be generated'.format(self.VID))

    # calculating the analytical stats: CTPS, CTRR.
    def analysis_stats(self):
        
        if self.CVNA < 2:
            self.CTPS = None
            self.CTRR = '.'
            return
        
        # logging.debug('CTRS score: {}'.format(self.CTRS))
        bins = [x[1] for x in g.ctps_cutoffs[:6]]
        bins = np.nextafter(bins, np.minimum(bins,0))
        ctps_index = np.digitize(self.CTRS, bins, right=True).tolist()
        self.CTPS = g.ctps_cutoffs[ctps_index][0]
        
        run_already = False
        cvcs_index = None
        for clinsig in self.CVCS.split(', '):
            if clinsig in [x[0] for x in g.ctps_cutoffs] and not run_already:
                run_already = True
                cvcs_index = [x[0] for x in g.ctps_cutoffs].index(clinsig)
            elif clinsig == 'Conflicting interpretations of pathogenicity':
                run_already = True
                cvcs_index = 3
            elif clinsig in [x[0] for x in g.ctps_cutoffs] and run_already:
                logging.warning('multiple CVCS statements for {}, only using '
                             'first one!'.format(self.VID))
        
        if cvcs_index is None:
            logging.warning('ClinVar significance for {} does not include B,B/LB'
                         ',LB,US,LP,LP/P,P'.format(self.VID))
            self.CTPS = None
            self.CTRR = '.'
            return
        
        self.CTRR = reclassification_tree(ctps_index, cvcs_index)
        
# test
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logging.debug('this is a module test')
    tree = ET.parse('../test/sample.xml')
    clinvarresult = tree.getroot()
    for var_index, variationreport in enumerate(clinvarresult):
        variant = VariationClass(variationreport)
        pprint.pprint(variant.__dict__)
    logging.debug('{} variant objects parsed in test'.format((var_index + 1)))
