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

from cogent.core.genetic_code import DEFAULT as standard_code
from cogent.db.ensembl import Genome,Species,host
import sqlalchemy as sql
import proBAM_biomart



#
# Retrieve correct Ensembl identifiers
#
def get_Ensembl_prefix(species):
    '''
    :param species: fulle species name
    :return: list of ENSEMBL transcript and protein tags
    '''
    ensembl_prefix=[]
    if species=='homo_sapiens':
        ensembl_prefix.append(['ENST','ENSP'])
    elif species=='mus_musculus':
        ensembl_prefix.append(['ENSMUST','ENSMUSP'])
    elif species=='drosophila_melanogaster':
        ensembl_prefix.append(['FBtr','FBpp'])
    elif species=='danio_rerio':
        ensembl_prefix.append(['ENSDART','ENSDARP'])
    elif species=='arabidopsis_thaliana':
        ensembl_prefix.append(['AT','AT'])
    else:
        raise ValueError('Species not recognized')
    return ensembl_prefix

#
# Retrieve data from BioMart and store in hashes
#
def prepareAnnotationENSEMBL(psm_protein_id,mode,database_v,species,three_frame_translation):
    '''
    :param psm_protein_id: list of protein IDs (untagged)
    :param mode: transcript or protein mode
    :param database_v: database version
    :param species: species name
    :return: dictionairy mapping proteins into ENSEMBL
    '''

    print('Commencing ENSEMBL data retrieval')
    # create connection to ensembl database
    if species=="arabidopsis_thaliana":
        Genome_species=Species.getCommonName(species.replace('_',' '))
        account=host.HostAccount(host="mysql-eg-publicsql.ebi.ac.uk",user="anonymous",passwd="",port=4157)
        ensembl=Genome(Species=Genome_species,account=account,Release=30)
    else:
        Genome_species=Species.getCommonName(species.replace('_',' '))
        ensembl=Genome(Species=Genome_species,Release=database_v,account=None)


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

    if mode=='protein':
        id=1
        query = sql.select(select_obj,from_obj=[from_obj],
                           whereclause = translation_table.c.stable_id.in_(psm_protein_id))

    elif mode=='transcript':
        id=0
        query = sql.select(select_obj,from_obj=[from_obj],
                           whereclause = transcript_table.c.stable_id.in_(psm_protein_id))
    psm_protein_id={}
    transcript_ids=[]
    for row in query.execute():
        #print row
        transcript_ids.append(row[2])
        psm_protein_id[row[id]]={'transcript_id':row[0],'translation_id':row[1],
                                 'transcript_seq':'','protein_seq':'',
                                 'chr':'','strand':'','5UTR_offset':row[3],'start_exon_rank':row[4]}
    return ensembl_construct_sequences(psm_protein_id,ensembl,transcript_ids,database_v,species,three_frame_translation)

#
# Retrieve DNA-sequence, Splice, construct protein sequence and store
#
def ensembl_construct_sequences(psm_hash,ensembl,transcript_ids,database_v,species,three_frame_translation):
    '''
    :param psm_hash: dictionair with protein / ensembl information ( see prepareAnnotationENSEMBL)
    :param ensembl:ensembl genome
    :param transcript_ids: list of transcrip ids (converted from protein IDs)
    :param database_v: database version
    :param species: species name
    :return: dictionairy mapping proteins into ENSEMBL
    '''
    print "Commencing transcript and protein sequence retrieval"

    no_protein_seq=[]
    biomart_key_hash={}
    stable_transcript_ids=[]

    for key in psm_hash.keys():
        biomart_key_hash[psm_hash[key]['transcript_id']]=key
        stable_transcript_ids.append(psm_hash[key]['transcript_id'])

    # Retrieve cds,chr,transcript_id and strand from biomart
    biomart_result=proBAM_biomart.retrieve_data_from_biomart(database_v,species,stable_transcript_ids,three_frame_translation)
    for row in biomart_result:
        row=row.split("\t")
        try:
            psm_hash[biomart_key_hash[row[1]]]['transcript_seq']=row[0]
            psm_hash[biomart_key_hash[row[1]]]['protein_seq']=standard_code.translate(row[0])
            #TODO what to do with "special" ensembl chromosomes: currently leave them out => bam conversion
            #TODO considers these psms unmapped
            #if "_" in row[2]:
            #    print row[1],row[2]
            psm_hash[biomart_key_hash[row[1]]]['chr']=row[2]
            psm_hash[biomart_key_hash[row[1]]]['strand']=row[3]
            del row
        except IndexError:
            pass
    del biomart_result

    # get exons directly from core database
    temp_exon_hash=get_ensembl_exons(ensembl,transcript_ids,psm_hash)
    exon_hash=temp_exon_hash[0]
    psm_hash=temp_exon_hash[1]
    del temp_exon_hash
    # retrieve protein sequences for transcript where the protein sequence could not be fetched automatically
    for key in no_protein_seq:
        psm_hash[key]['transcript_seq']=retrieve_protein_seq(psm_hash[key]['transcript_seq'],
                                                                   exon_hash[psm_hash[key]['transcript_id']],
                                                                   psm_hash[key]['5UTR_offset'],
                                                                   psm_hash[key]['start_exon_rank'])
        #translate till stop codon
        psm_hash[key]['protein_seq']=standard_code.translate(psm_hash[key]
                                                                   ['transcript_seq']).partition('*')[0]
    return [psm_hash,exon_hash]

#
# Retrieve protein coding sequence from transcript
#

def retrieve_protein_seq(transcript_seq,exons,offset,start_exon):
    '''
    :param transcript_seq: transcript sequence
    :param exons: exons from transcript sequence
    :param offset: offset from transcript start site till CDS start site
    :param start_exon: first exon containing CDS
    :return: CDS sequence
    '''
    count=0
    for exon in exons:
        exon=[int(numeric_string) for numeric_string in exon]
        if exon[2]<start_exon:
            count+=exon[1]-exon[0]
        elif exon[2]==start_exon:
            count+=offset
            break
    return transcript_seq[(count-1):len(transcript_seq)]

def get_ensembl_exons(ensembl,transcript_ids,psm_hash):
    '''
    :param ensembl: ENSEMBL genome
    :param transcript_ids: list of transcript ids
    :param psm_hash: dictionairy mapping proteins as transcript on ensembl
    :return: exon hash
    '''
    print "Commencing exon retrieval"
    exon_hash={}
    #connect and join tables in order to retrieve exon information
    transcript_table=ensembl.CoreDb.getTable('transcript')
    exon_transcript_table=ensembl.CoreDb.getTable('exon_transcript')
    exon_table=ensembl.CoreDb.getTable('exon')
    select_obj=[exon_table.c.seq_region_start,
                exon_table.c.seq_region_end,
                exon_transcript_table.c.rank,
                transcript_table.c.stable_id,
                exon_transcript_table.c.exon_id
                ]
    from_obj=exon_transcript_table.join(exon_table,exon_table.c.exon_id==exon_transcript_table.c.exon_id) \
                                  .join(transcript_table,transcript_table.c.transcript_id==
                                        exon_transcript_table.c.transcript_id)
    query = sql.select(select_obj,from_obj=[from_obj],
                           whereclause = exon_transcript_table.c.transcript_id.in_(transcript_ids))

    #store exon information in hash, return hash
    for row in query.execute():
        if row[3] not in exon_hash.keys():
            exon_hash[row[3]]=[]
        if row[3] in psm_hash.keys():
            if psm_hash[row[3]]['start_exon_rank']==row[4]:
                psm_hash[row[3]]['start_exon_rank']=row[2]
        exon_hash[row[3]].append([str(row[0]),str(row[1]),str(row[2])])
    return [exon_hash,psm_hash]

#
# Create SQ header for SAM file (ENSEMBL) get chr location and coordinates and store them in an array
#
def create_SQ_header(database_v,species):
    '''
    :param database_v: database version
    :return: list of chromosomes with their size (from ENSEMBL)
    '''
    SQ=[]
    # create connection to ensembm database
        # create connection to ensembl database
    if species=="arabidopsis_thaliana":
        Genome_species=Species.getCommonName(species.replace('_',' '))
        account=host.HostAccount(host="mysql-eg-publicsql.ebi.ac.uk",user="anonymous",passwd="",port=4157)
        ensembl=Genome(Species=Genome_species,account=account,Release=30)
    else:
        Genome_species=Species.getCommonName(species.replace('_',' '))
        ensembl=Genome(Species=Genome_species,Release=database_v,account=None)

    # convert IDs
    coord_table=ensembl.CoreDb.getTable('coord_system')
    seq_region_table=ensembl.CoreDb.getTable('seq_region')
    select_obj=[seq_region_table.c.name,
                seq_region_table.c.length,
                coord_table.c.version,
                ]
    from_obj=seq_region_table.join(coord_table,coord_table.c.coord_system_id==seq_region_table.c.coord_system_id)

    query = sql.select(select_obj,from_obj=[from_obj],
                           whereclause = coord_table.c.rank==1)

    for row in query.execute():
        if '_' not in row[0]:
            SQ_string= "@SQ\tSN:chr"+str(row[0])+"\tLN:"+str(row[1])+"\tAS:"+str(row[2])+"\tSP:"+str(species)
            SQ.append(SQ_string)
    return SQ