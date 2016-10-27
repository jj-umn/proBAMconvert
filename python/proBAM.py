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

from __future__ import division

__author__ = 'Volodimir Olexiouk'

import time
from itertools import imap
import operator
import re
import pysam
import pysam.ctabixproxies
import proBAM_biomart
import proBAM_input
import proBAM_ENSEMBL
import sys
import getopt
import os
import proBAM_IDparser
import proBAM_pepxml
from cogent.core.genetic_code import DEFAULT as standard_code
from proBAM_coref import *

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

######################################
###         DEPENDENCIES           ###
######################################
'''
cogent
pysam
mySQLdb
SQLAlchemy
lxml
numpy
matplotlib
pyteomics
'''
'''
#
# Command line variable input
#
def get_input_variables():

    ###############################
    # GLOBAL VARIABLE DECLARATION #
    ###############################

    global name
    global allowed_mismatches
    global database_v
    global database
    global species
    global psm_file
    global directory
    global decoy_annotation
    global version
    global sorting_order
    global map_decoy
    global rm_duplicates
    global three_frame_translation

    ######################################
    ### VARIABLE AND INPUT DECLARATION ###
    ######################################

    name=""
    mismatches=0
    version=""
    database=""
    species=""
    psm_file=""
    directory=""
    decoy_annotation=['REV_','DECOY_','_REVERSED']
    version='1.0'
    # can be unknown,unsorted, queryname or coordinate, can be specified by user
    sorting_order='unknown'
    map_decoy=""
    rm_duplicates=""
    three_frame_translation


    #
    # Read command line args
    # Check for valid arguments
    #

    try:
        myopts, args = getopt.getopt(sys.argv[1:],"n:m:v:d:s:f:c:r:a:t",["name=","mismatches=","version=","database=","species=","file=","directory=","rm_duplicates=","map_decoy=","tri_frame_translation="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit()

    ###############################
    # o == option
    # a == argument passed to the o
    ###############################

    #
    # Catch arguments
    #

    for o, a in myopts:
        if o in ('-n','--name'):
            result_db=a
        if o in ('-m','--mismatches'):
            mismatches=a
        if o in ('-v','--version'):
            version=a
        if o in ('-d','--database'):
            database=a
        if o in ('-s','--species'):
            species=a
        if o in ('-f','--file'):
            psm_file=a
        if o in ('-c','--directory'):
            directory=a
        if o in ('-r','--rm_duplicates'):
            rm_duplicates=a
        if o in ('-a','--map_decoy'):
            map_decoy=a
        if o in ('-t','--three_frame_translation'):
            map_decoy=a

    #
    # Check for correct argument, output argument and parse
    #

    if(database==''):
        print("Error: do not forget to pass the database")
        sys.exit()
    database=database.upper()
    if database !="ENSEMBL":
        print 'Error: unsupported database \n' \
              'currently supported databases: ENSEMBL'
        sys.exit
    if(version == ''):
        print("Error: do not forget to pass the database version !")
        sys.exit()
    species=species.lower()
    proBAM_biomart._get_ensembl_dataset_(species)
    if(species==''):
        print("Error: do not forget to pass the species argument")
        sys.exit()
    if(file == ''):
        print("Error: do not forget to pass psm file (pepxml, mzid or mztab) !")
        sys.exit()
    if(directory == ''):
       directory = os.getcwd()
       os.chdir(directory)
    if(directory !=''):
        os.chdir(directory)
    if rm_duplicates!="Y":
        rm_duplicats="N"
    if map_decoy!="Y":
        map_decoy="N"
    if three_frame_translation!="Y":
        three_frame_translation="N"

    allowed_mismatches=mismatches
    database_v=version

    command_line= python proBAM.py --name=,"mismatches=","version=","database=","species=","file=","directory=","rm_duplicates=","map_decoy=","tri_frame_translation="


    # ouput variables
    print(  "psm file:                                      " + str(file) +"\n"+
            "directory:                                     " + str(directory) +"\n"+
            "database                                       " + str(database) +"\n"+
            "database version:                              " + str(database_v) +"\n"+
            "species:                                       " + str(species) +"\n"+
            "allowed mismatches:                            " + str(allowed_mismatches)+"\n"+
            "three_frame_translation:                         " + str(three_frame_translation))
'''


###############################
# NON COMMAND LINE ARGUMENTS  #
# FOR TESTING PURPOSES        #
###############################

directory="/home/vladie/Desktop/proBAMconvert/output/"
psm_file="/home/vladie/Desktop/proBAMconvert/PeptideShaker.mzid"
species="homo_sapiens"
database='ENSEMBL'
database_v=85
# TODO Let users specify used the decoy annotation
decoy_annotation=['REV_','DECOY_','_REVERSED']
allowed_mismatches=1
version='1.0'
# can be unknown,unsorted, queryname or coordinate, can be specified by user
sorting_order='unknown'
name='test'
three_frame_translation='N'
map_decoy="N"
rm_duplicates="N"

command_line= "python proBAM.py --name "+str(name)+" --mismatches "+str(allowed_mismatches)+" --version "+str(database_v)\
              +" --database "+str(database)+" --species "+str(species)+" --file "+str(psm_file)+\
              " --directory "+str(directory)+" --rm_duplicates "+str(rm_duplicates)+\
              " --map_decoy "+str(map_decoy)+" --tri_frame_translation "+str(three_frame_translation)

# ouput variables
print(  "psm file:                                      " + str(psm_file) +"\n"+
        "directory:                                     " + str(directory) +"\n"+
        "database                                       " + str(database) +"\n"+
        "database version:                              " + str(database_v) +"\n"+
        "species:                                       " + str(species) +"\n"+
        "allowed mismatches:                            " + str(allowed_mismatches)+"\n"+
        "three_frame_translation:                       " + str(three_frame_translation)+"\n"+
        "map decoys:                                    " + str(map_decoy)+"\n"+
        "remove duplicates:                             " + str(rm_duplicates))


#######################
### GETTERS/SETTERS ###
#######################

#
# Function to write temp_result to file
#
def write_psm(temp_result,file):
    '''
    :param temp_result: sam results
    :param file: sam file
    :return: write to sam file IO
    '''
    for tab in temp_result:
        file.write(str(tab)+'\t')
    file.write('\n')

#
# Open file to write to
#
def open_sam_file(directory,name):
    '''
    :param directory: working directory
    :param name: output file name
    :return: return output file IO
    '''
    file=open(directory+name+'.sam',"w")
    return file


######################
### MAIN FUNCTIONS ###
######################

#
# Convert PSM to SAM
#
def PSM2SAM(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,map_decoy,rm_duplicates,
            three_frame_translation,psm_file):
    '''
    :param psm_hash: dictionairy of psm files
    :param transcript_hash: dictionairy of transcripts
    :param exon_hash: dictionairy of exons
    :param decoy_annotation: decoy annotation list
    :param allowed_mismatches: number of allowed mismatches
    :param file: sam file
    :return: sam file IO
    '''
    # psm_hash.reset()
    if rm_duplicates == "Y":
        dup = {}
    enzyme=get_enzyme(psm_file)
    enzyme_specificity=get_enzyme_specificity(psm_file)
    for psm in psm_hash:
        # convert unmapped PSMs with their own converter
        if 'search_hit' not in psm.keys():
            nohit_PSM_to_SAM(psm,file,enzyme,enzyme_specificity)
        else:
            for row in psm['search_hit']:
                for p in range(0,len(row['proteins'])):
                    decoy=0
                    # convert decoys with decoy-specific convertor
                    for d in decoy_annotation:
                        if d in row['proteins'][p]['protein'].upper():
                            decoy=1
                            key=row['proteins'][p]['protein'].upper().split(d)[1]
                            if key in transcript_hash.keys():
                                temp_result= decoy_PSM_to_SAM(psm,row,key,transcript_hash,exon_hash,allowed_mismatches,
                                                 map_decoy,enzyme,enzyme_specificity)
                                if rm_duplicates=="Y":
                                    dup_key= str(temp_result[9])+"_"+str(temp_result[2])+"_"+str(temp_result[3])
                                    if dup_key not in dup.keys():
                                        dup[dup_key]=1
                                        write_psm(temp_result,file)
                                else:
                                    write_psm(temp_result,file)
                            else:
                                temp_result=unannotated_PSM_to_SAM(psm,row,decoy,key,enzyme,enzyme_specificity)
                                if rm_duplicates=="Y":
                                    dup_key= str(temp_result[9])+"_"+str(temp_result[2])+"_"+str(temp_result[3])
                                    if dup_key not in dup.keys():
                                        dup[dup_key]=1
                                        write_psm(temp_result,file)
                                else:
                                    write_psm(temp_result,file)

                    if decoy==0:
                        key=row['proteins'][p]['protein']
                        # Filter out PSM where transcript sequences were not found/ non-existent
                        if key not in transcript_hash.keys():
                            write_psm(unannotated_PSM_to_SAM(psm,row,decoy,key,enzyme,enzyme_specificity),file)
                        else:
                            if three_frame_translation=="Y":
                                protein_hit=map_peptide_to_protein_3frame(row['peptide'],
                                                                          transcript_hash[key]['transcript_seq'],
                                                                          allowed_mismatches,
                                                                          transcript_hash[key]['strand'])[0]
                                pre_post_aa=map_peptide_to_protein_3frame(row['peptide'],
                                                                          transcript_hash[key]['transcript_seq'],
                                                                          allowed_mismatches,
                                                                          transcript_hash[key]['strand'])[1]
                            else:
                                protein_hit=map_peptide_to_protein(row['peptide'],transcript_hash[key]['protein_seq']
                                                                   ,allowed_mismatches)[0]
                                pre_post_aa=map_peptide_to_protein(row['peptide'],transcript_hash[key]['protein_seq']
                                                                   ,allowed_mismatches)[1]
                            if len(protein_hit)==0:
                                write_psm(unannotated_PSM_to_SAM(psm,row,decoy,key,enzyme,enzyme_specificity),file)
                            else:
                                # map peptide on protein and retrieve hit position, iterate over all hits
                                for phit in protein_hit:
                                    temp_result=[None]*32
                                    #
                                    # Mandatory columns adapted from SAM/BAM format
                                    #
                                    #QNAME
                                    temp_result[0]=psm['spectrum']
                                    #FLAG
                                    temp_result[1]=calculate_FLAG(transcript_hash[key]['strand'],row['hit_rank'],
                                                                  decoy)
                                    #RNAME
                                    temp_result[2]='chr'+str(transcript_hash[key]['chr'])
                                    #POS
                                    temp_result[3]=calculate_genome_position(phit[0],
                                                                             transcript_hash[key]['strand'],
                                                                             transcript_hash[key]['5UTR_offset'],
                                                                             transcript_hash[key]['start_exon_rank'],
                                                                             row['peptide'],
                                                                             exon_hash[transcript_hash[key]['transcript_id']],
                                                                             transcript_hash[key]['chr'],
                                                                             three_frame_translation)
                                    #MAPQ
                                    temp_result[4]=255
                                    #CIGAR
                                    temp_result[5]=compute_cigar(temp_result[3],
                                                                 exon_hash[transcript_hash[key]['transcript_id']],
                                                                 transcript_hash[key]['strand'],row['peptide'])[0]
                                    #RNEXT
                                    temp_result[6]='*'
                                    #PNEXT
                                    temp_result[7]=0
                                    #TLEN
                                    temp_result[8]=0
                                    #SEQ
                                    if three_frame_translation=='Y':
                                        phit_loc=phit[0]
                                    else:
                                        phit_loc=phit[0]*3
                                    if int(transcript_hash[key]['strand'])==1:
                                        temp_result[9]=str(transcript_hash[key]['transcript_seq']\
                                                   [phit_loc:(phit_loc+(len(row['peptide'])*3))])
                                    else:
                                        temp_result[9]=reverse_complement(str(transcript_hash[key]['transcript_seq']\
                                                   [phit_loc:(phit_loc+(len(row['peptide'])*3))]))
                                    #QUAL
                                    temp_result[10]='*'
                                    #
                                    #Mandatory proteomics specific columns added to the proBam format
                                    #
                                    #NH: number of genomic location the peptide mapping to
                                    temp_result[11]='NH:i:'+str(len(row['proteins'])+len(phit)-1)
                                    #XO: uniqueness of peptide mapping
                                    #todo figure this one out
                                    temp_result[12]='XO:Z:*'
                                    #XL: number of peptides the spectrum mapping to
                                    temp_result[13]='XL:i:'+str(len(psm['search_hit']))
                                    #XP; peptide sequence
                                    temp_result[14]='XP:Z:'+row['modified_peptide']
                                    #YP: protein accession ID from the original search
                                    temp_result[15]='YP:Z:'+str(key)
                                    #XF: reading frame of the peptide
                                    temp_result[16]='XF:Z:'+compute_cigar(temp_result[3],
                                                                 exon_hash[transcript_hash[key]['transcript_id']],
                                                                 transcript_hash[key]['strand'],row['peptide'])[1]
                                    #XI: peptide intensity
                                    temp_result[17]="XI:f:*"
                                    #XB: Mass error (experimental - calculated)
                                    temp_result[18]="XB:f:"+str(row['massdiff'])
                                    #XR: reference peptide sequence
                                    temp_result[19]='XR:Z:'+row['peptide']
                                    #YB: preceding 2AA
                                    temp_result[20]="YB:Z:"+str(pre_post_aa[0])
                                    #YA: following 2AA:
                                    temp_result[21]="YA:Z:"+str(pre_post_aa[1])
                                    #XS: PSM score
                                    temp_result[22]=create_XS(row['search_score'])
                                    #XQ: PSM-Qvalue
                                    temp_result[23]='XQ:f:'+'NA'
                                    #XC: Peptide charge
                                    temp_result[24]='XC:i:'+str(psm['assumed_charge'])
                                    #XA: Whether the peptide is well annotated
                                    temp_result[25]=create_XA(phit[1])
                                    #XM: Modification
                                    temp_result[26]=create_XM(row['modifications'])
                                    #XN: number of mis-cleavages
                                    if 'num_missed_cleaveges' in row.keys():
                                        temp_result[27]='XN:i:'+str(row['num_missed_cleavages'])
                                    else:
                                        temp_result[27]='XN:i:0'
                                    #XT: non/semi/tryptic
                                    temp_result[28]="XT:i:"+str(enzyme_specificity)
                                    #XE: enzyme used
                                    temp_result[29]="XE:i:"+str(enzyme)
                                    #XG: Petide type
                                    temp_result[30]=create_XG(phit[1])
                                    #XU: petide URL
                                    temp_result[31]="XU:Z:*"
                                    # remove duplicates if rm_duplicates=Y
                                    if rm_duplicates=="Y":
                                        dup_key= str(temp_result[9])+"_"+\
                                                          str(temp_result[2])+"_"+str(temp_result[3])
                                        if dup_key not in dup.keys():
                                            dup[dup_key]=1
                                            write_psm(temp_result,file)
                                    else:
                                        write_psm(temp_result,file)
    file.close()

#
# Function to convert unmapped PSM to SAM
#

def nohit_PSM_to_SAM(psm,file,enzyme,enzyme_specificity):
    '''
    :param psm: psm file
    :param file: output file
    :return: converted psm to SAM if no hits found
    '''
    temp_result=[None]*32
    #
    # Mandatory columns adapted from SAM/BAM format
    #
    #QNAME
    temp_result[0]=psm['spectrum']
    #FLAG
    temp_result[1]='4'
    #RNAME
    temp_result[2]='*'
    #POS
    temp_result[3]=0
    #MAPQ
    temp_result[4]=255
    #CIGAR
    temp_result[5]='*'
    #RNEXT
    temp_result[6]='*'
    #PNEXT
    temp_result[7]=0
    #TLEN
    temp_result[8]=0
    #SEQ
    temp_result[9]='*'
    #QUAL
    temp_result[10]='*'
    #
    #Mandatory proteomics specific columns added to the proBam format
    #
    #NH: number of genomic location the peptide mapping to
    temp_result[11]='NH:i:0'
    #XO uniqueness of peptide mapping
    temp_result[12]='XO:Z:*'
    #XL: number of peptides the spectrum mapping to
    temp_result[13]='XL:i:0'
    #XP; peptide sequence
    temp_result[14]='XP:Z:NA'
    #YP: Protein accession ID from original search
    temp_result[15]='YP:Z:*'
    #XF : reading frame of the peptide
    temp_result[16]='XF:i:*'
    #XI : Peptide intensity
    temp_result[17]='XI:Z:*'
    #XB: mass error
    temp_result[18]='XB:f:*'
    #XR: reference peptide sequence
    temp_result[19]='XR:Z:NA'
    #YB: 2AA before
    temp_result[20]='YB:Z:*'
    #YA: 2AA after
    temp_result[21]='YA:Z:*'
    #XS: PSM score
    temp_result[22]='XS:f:0'
    #XQ: PSM-Qvalue
    temp_result[23]='XQ:f:0'
    #XC: Peptide charge
    temp_result[24]='XC:i:0'
    #XA: Whether the peptide is well annotated
    temp_result[25]='XA:i:0'
    #XM: Modification
    temp_result[26]='XM:Z:NA'
    #XN: number of mis-cleavages
    temp_result[27]='XN:i:0'
    #XT: non/semi/tryptic
    temp_result[28]="XT:i:"+str(enzyme_specificity)
    #XE: enzyme used
    temp_result[29]="XE:i:"+str(enzyme)
    #XG: Petide type
    temp_result[30]="XG:Z:U"
    #XU: URL
    temp_result[31]="XU:Z:*"
    write_psm(temp_result,file)


#
# Function to convert unannotated PSMs to SAM
#
def unannotated_PSM_to_SAM(psm,row,decoy,key,enzyme,enzyme_specificity):
    '''
    :param psm: psm dictionairy
    :param row: unnanotated PSM row
    :param decoy: decoy boolean
    :param file: output file
    :return: sam of unnanotated PSM
    '''
    decoy=int(decoy)
    temp_result=[None]*32
    #
    # Mandatory columns adapted from SAM/BAM format
    #
    #QNAME
    temp_result[0]=psm['spectrum']
    #FLAG
    if decoy==1:
        temp_result[1]='4'
    else:
        temp_result[1]='1028'
    #RNAME
    temp_result[2]='*'
    #POS
    temp_result[3]=0
    #MAPQ
    temp_result[4]=255
    #CIGAR
    temp_result[5]='*'
    #RNEXT
    temp_result[6]='*'
    #PNEXT
    temp_result[7]=0
    #TLEN
    temp_result[8]=0
    #SEQ
    temp_result[9]='*'
    #QUAL
    temp_result[10]='*'
    #
    #Mandatory proteomics specific columns added to the proBam format
    #
    #NH: number of genomic location the peptide mapping to
    temp_result[11]='NH:i:'+str(len(row['proteins']))
    #XO: uniqueness
    temp_result[12]='XO:Z:*'
    #XL: number of peptides the spectrum mapping to
    temp_result[13]='XL:i:'+str(len(psm['search_hit']))
    #XP; peptide sequence
    temp_result[14]='XP:Z:'+row['modified_peptide']
    #YP: protein accession id
    temp_result[15]="YP:Z:"+str(key)
    #XF: Reading frame of the peptide
    temp_result[16]="XF:Z:*"
    #XI: Peptide intensity
    temp_result[17]="XI:f:*"
    #XB: Mass error
    temp_result[18]="XB:f:"+str(row['massdiff'])
    #XR: reference peptide sequence
    temp_result[19]='XR:Z:'+row['peptide']
    #YB: 2 AA before
    temp_result[20]='YB:Z:*'
    #YA: 2 AA after
    temp_result[21]='YA:Z:*'
    #XS: PSM score
    temp_result[22]=create_XS(row['search_score'])
    #XQ: PSM-Qvalue
    temp_result[23]='XQ:f:0'
    #XC: Peptide charge
    temp_result[24]='XC:i:'+str(psm['assumed_charge'])
    #XA: Whether the peptide is well annotated
    temp_result[25]='XA:i:2'
    #XM: Modification
    temp_result[26]=create_XM(row['modifications'])
    #XN: number of mis-cleavages
    if 'num_missed_cleavages' in row.keys():
        temp_result[27]='XN:i:'+str(row['num_missed_cleavages'])
    else:
        temp_result[27]='XN:i:0'
    #XT: non/semi/tryptic
    temp_result[28]="XT:i:"+str(enzyme_specificity)
    #XE: enzyme
    temp_result[29]="XE:i:"+str(enzyme)
    #XG: Petide type
    if decoy==1:
        temp_result[30]="XG:Z:U"
    else:
        temp_result[30]="XG:Z:D"
    #XU
    temp_result[31]="XU:Z:*"
    return temp_result

#
# Function to convert decoy PSM to SAM format
#

def decoy_PSM_to_SAM(psm,row,key,transcript_hash,exon_hash,allowed_mismatches,map_decoy,enzyme,enzyme_specificity):
    '''
    :param psm: psm dictionairy
    :param row: row where decoy found
    :param key: psm key
    :param transcript_hash: transcript dictionairy
    :param exon_hash: exon dictionairy
    :param allowed_mismatches: number of allowed mismatches
    :param file: output file
    :return: SAM of decoy PSM
    '''
    temp_result=[None]*23
    #map decoy to genome if map_decoy=Y
    if map_decoy=="Y":
        protein_hit=map_peptide_to_protein(row['peptide'][::-1],transcript_hash[key]['protein_seq'],allowed_mismatches)[0]
        pre_post_aa=map_peptide_to_protein(row['peptide'][::-1],transcript_hash[key]['protein_seq'],allowed_mismatches)[1]
    else:
        protein_hit=[]
    if len(protein_hit)==0:
         return unannotated_PSM_to_SAM(psm,row,1,key,enzyme,enzyme_specificity)
    else:
        # map peptide on protein and retrieve hit position, iterate over all hits
        for phit in protein_hit:
            temp_result=[None]*32
            #
            # Mandatory columns adapted from SAM/BAM format
            #
            #QNAME
            temp_result[0]=psm['spectrum']
            #FLAG
            temp_result[1]=calculate_FLAG(transcript_hash[key]['strand'],row['hit_rank'],
                                          1)
            #RNAME
            temp_result[2]='chr'+str(transcript_hash[key]['chr'])
            #POS
            temp_result[3]=calculate_genome_position(phit[0],
                                                     transcript_hash[key]['strand'],
                                                     transcript_hash[key]['5UTR_offset'],
                                                     transcript_hash[key]['start_exon_rank'],
                                                     row['peptide'][::-1],
                                                     exon_hash[transcript_hash[key]['transcript_id']],
                                                     transcript_hash[key]['chr'])
            #MAPQ
            temp_result[4]=255
            #CIGAR
            temp_result[5]=compute_cigar(temp_result[3],
                                         exon_hash[transcript_hash[key]['transcript_id']],
                                         transcript_hash[key]['strand'],row['peptide'])
            #RNEXT
            temp_result[6]='*'
            #PNEXT
            temp_result[7]=0
            #TLEN
            temp_result[8]=0
            #SEQ
            if int(transcript_hash[key]['strand'])==1:
                temp_result[9]=str(transcript_hash[key]['transcript_seq']\
                           [(phit[0]*3):((phit[0]*3)+(len(row['peptide'])*3))])
            else:
                temp_result[9]=reverse_complement(str(transcript_hash[key]['transcript_seq']\
                           [(phit[0]*3):((phit[0]*3)+(len(row['peptide'])*3))]))
            #QUAL
            temp_result[10]='*'
            #
            #Mandatory proteomics specific columns added to the proBam format
            #
            #NH: number of genomic location the peptide mapping to
            temp_result[11]='NH:i:'+str(len(row['proteins'])+len(phit)-1)
            # todo figure this one out
            temp_result[12] = 'XO:z:*'
            # XL: number of peptides the spectrum mapping to
            temp_result[13] = 'XL:i:' + str(len(psm['search_hit']))
            # XP; peptide sequence
            temp_result[14] = 'XP:Z:' + row['modified_peptide']
            # YP: protein accession ID from the original search
            temp_result[15] = 'YP:Z:' + str(key)
            # XF: reading frame of the peptide
            temp_result[16] = 'XF:Z:' + compute_cigar(temp_result[3],
                                                      exon_hash[transcript_hash[key]['transcript_id']],
                                                      transcript_hash[key]['strand'], row['peptide'])[1]
            # XI: peptide intensity
            temp_result[17] = "XI:f:*"
            # XB: Mass error (experimental - calculated)
            temp_result[18] = "XB:f:" + str(row['massdiff'])
            # XR: reference peptide sequence
            temp_result[19] = 'XR:Z:' + row['peptide']
            # YB: preceding 2AA
            temp_result[20] = "YB:Z:" + str(pre_post_aa[0])
            # YA: following 2AA:
            temp_result[21] = "YA:Z:" + str(pre_post_aa[1])
            #XS: PSM score
            temp_result[22]=create_XS(row['search_score'])
            #XQ: PSM-Qvalue
            temp_result[23]='XQ:f:'+'NA'
            #XC: Peptide charge
            temp_result[24]='XC:i:'+str(psm['assumed_charge'])
            #XA: Whether the peptide is well annotated
            temp_result[25]=create_XA(phit[1])
            #XM: Modification
            temp_result[26]=create_XM(row['modifications'])
            #XN: number of mis-cleavages
            if 'num_missed_cleaveges' in row.keys():
                temp_result[27]='XN:i:'+str(row['num_missed_cleavages'])
            else:
                temp_result[27]='XN:i:0'
            #XT: non/semi/tryptic
            temp_result[28]="XT:i:"+str(enzyme_specificity)
            #XE enzyme
            temp_result[29]="XE:i"+str(enzyme)
            #XG: Petide type
            temp_result[30]='XG:Z:D'
            #XU= url
            temp_result[31]="XU:Z:*"
            return temp_result

#
# Create SAM header
#
def create_SAM_header(file,version,database,sorting_order,database_v,species,command_line,psm_file):
    '''
    :param file: output file
    :param version: proBAMconvert version
    :param database: database name
    :param sorting_order: SAM sorting order
    :param database_v: database version
    :return:
    '''
    header=[]
    header.append('@HD\tVN:'+version+' SO:'+sorting_order)
    if database=="ENSEMBL":
        SQ=proBAM_ENSEMBL.create_SQ_header(database_v,species)
    for row in SQ:
        header.append(row)
    header.append('@PG\tID:proBamPy\tVN:1.0\tCL:'+str(command_line))
    header.append('@GA\tAS:'+str(database)+'\tVN:'+str(database_v))
    # get comments and append comments to file
    comments=extract_comments(psm_file)
    for comment in comments:
        comment=str(comment).rstrip()
        comment=re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', comment, flags=re.M)
        header.append('@CO\t'+str(comment))
    for row in header:
        file.write(row+'\n')


#
# Function to convert SAM to BAM
#
def sam_2_bam(directory,name):
    '''
    :param directory: working directory
    :param name: file name
    :return: bam file
    '''
    infile = pysam.AlignmentFile((directory+name+'.sam'), "r")
    outfile = pysam.AlignmentFile((directory+name+'.bam'), "wb", template=infile)
    for s in infile:
        outfile.write(s)
    pysam.sort((directory+name+'.bam'),"-o",(directory+name+'.sorted.bam'))
    pysam.index(directory+name+'.sorted.bam')



####################
### MAIN PROGRAM ###
####################

if __name__=='__main__':
    #get_input_variables()
    start_time = time.time()                                 # start timing function
    file=open_sam_file(directory,name)

    # hash PSM_DATA and define variables
    psm_hash=proBAM_input.get_PSM_hash(psm_file,decoy_annotation)
    parse_results=proBAM_IDparser.parseID(psm_hash,species,database,decoy_annotation,database_v,three_frame_translation)

    annotation=parse_results[1]
    psm_hash=parse_results[0]
    transcript_hash=annotation[0]
    exon_hash=annotation[1]

    # convert to SAM
    create_SAM_header(file,version,database,sorting_order,database_v,species,command_line,psm_file)
    PSM2SAM(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,map_decoy,rm_duplicates,
            three_frame_translation,psm_file)
    sam_2_bam(directory,name)


    print("proBAM conversion succesful")
    print("%f seconds" % (time.time() - start_time))         # output script run time
