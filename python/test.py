__author__ = 'vladie'

from bioservices import *


biomart = BioMart(host="jul2015.archive.ensembl.org")
test=biomart.datasets("ENSEMBL_MART_ENSEMBL")
for dataset in test:
    print dataset