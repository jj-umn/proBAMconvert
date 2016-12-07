##############################################################################
# Copyright 2016  Olexiouk Volodimir,Menschaert Gerben                       #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# You may obtain a copy of the License at                                    #
#                                                                            #
#  http://www.apache.org/licenses/LICENSE-2.0                                #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
##############################################################################

__author__ = 'vladie'

from bioservices import *
from cogent.core.genetic_code import DEFAULT as standard_code
from cogent.db.ensembl import Genome,Species
import sqlalchemy as sql
import proBAM_biomart
import cogent.db.ensembl.host as test

'''
biomart=BioMart(host="plants.ensembl.org")
test=biomart.datasets('plants_mart')

for i in test:
    print i

#print Species
Genome_species=Species.getCommonName("arabidopsis_thaliana".replace('_',' '))
#print Genome_species
account=test.HostAccount(host="mysql-eg-publicsql.ebi.ac.uk",user="anonymous",passwd="",port=4157)
ensembl=Genome(Species="A.thaliana",account=account,release=30)




db=test.get_db_name(account=account,DEBUG=True,species="A.thaliana",db_type="core",release=30)
for i in db: print i

# convert IDs
translation_table=ensembl.CoreDb.getTable('translation')
transcript_table=ensembl.CoreDb.getTable('transcript')
select_obj=[transcript_table.c.stable_id,
            translation_table.c.stable_id,
            transcript_table.c.transcript_id,
            translation_table.c.seq_start,
            translation_table.c.start_exon_id,
            ]
from_obj=translation_table.join(transcript_table,transcript_table.c.transcript_id==translation_table.c.transcript_id)
query = sql.select(select_obj,from_obj=[from_obj])

psm_protein_id={}
transcript_ids=[]

for row in query.execute():
    print row



biomart=BioMart(host="plants.ensembl.org")
test2=biomart.attributes('athaliana_eg_gene')
print biomart.databases
biomart.configuration("athaliana_eg_gene")

biomart.add_dataset_to_xml('athaliana_eg_gene')
biomart.add_attribute_to_xml("ensembl_transcript_id")
biomart.add_attribute_to_xml("transcript_start")
biomart.add_attribute_to_xml("uniprot_swissprot_accession")
biomart.add_attribute_to_xml("transcript_end")
xml_query=biomart.get_xml()
xml_query=xml_query.replace('virtualSchemaName = "default"','virtualSchemaName = "plants_mart_30"')

temp_result=biomart.query(xml_query).split("\n")
result=[]
for row in temp_result:
    items=row.split("\t")
    print row
    if len(items)==4:
        length=int(items[3])-int(items[1])+1
        result.append(items[0]+" "+str(length)+" "+items[3])


from proBAM import map_peptide_to_protein_3frame

map_peptide_to_protein_3frame(' I L L Q E I G T','ATTTAGGACCAATAAGTCTTAATTGGTTTGAAGAACTTTCTTCAGAAGCTCCACCCTATAATTCTGAACCTGCAGAAGAATC' \
                               'TGAACATAAAAACAACAATTACGAACCAAACCTATTTAAAACTCCACAAAGGAAACCATCTTATAATCAGCTGGCTTCAACT' \
                               'CCAATAATATTCAAAGAGCAAGGGCTGACTCTGCCGCTGTACCAATCTCCTGTAAAAGAATTAGATAAATTCAAATTAGAC' \
                               'TTAG ',0,-1)
'''
#import re
#([OPQ][0-9][A-Z0-9]{3}[0-9])|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]{1,2}
#test=re.compile("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]")
#print re.findall(test,"A0A000")

import urllib,urllib2

def get_updated_entry_name(name):
    url = 'http://www.uniprot.org/uniprot/'+name+'.tab'

    params = {
    'columns':'id',
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "volodimir.olexiouk@ugent.be" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read(200000)
    return page.split("\n")[1]