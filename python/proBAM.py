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
from cogent.core.genetic_code import DEFAULT as standard_code

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


    # ouput variables
    print(  "psm file:                                      " + str(file) +"\n"+
            "directory:                                     " + str(directory) +"\n"+
            "database                                       " + str(database) +"\n"+
            "database version:                              " + str(database_v) +"\n"+
            "species:                                       " + str(species) +"\n"+
            "allowed mismatches:                            " + str(allowed_mismatches)+"\n"+
            "three_frame_translation:                         " + str(three_frame_translation))


###############################
# NON COMMAND LINE ARGUMENTS  #
# FOR TESTING PURPOSES        #
###############################
'''
directory="/home/vladie/Desktop/proBAM_testing/examples/"
psm_file="/home/vladie/Desktop/proBAM_testing/examples/Maria_Nano-particles_CSNP2_4of10ul_02.pep.xml"
species="homo_sapiens"
database='ENSEMBL'
database_v=82
# TODO Let users specify used the decoy annotation
decoy_annotation=['REV_','DECOY_','_REVERSED']
allowed_mismatches=0
version='1.0'
# can be unknown,unsorted, queryname or coordinate, can be specified by user
sorting_order='unknown'
name='test'
'''

#######################
### GETTERS/SETTERS ###
#######################



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
            three_frame_translation):
    '''
    :param psm_hash: dictionairy of psm files
    :param transcript_hash: dictionairy of transcripts
    :param exon_hash: dictionairy of exons
    :param decoy_annotation: decoy annotation list
    :param allowed_mismatches: number of allowed mismatches
    :param file: sam file
    :return: sam file IO
    '''
    print 'Converting pepXML to SAM format'
    # psm_hash.reset()
    for psm in psm_hash:

        # convert unmapped PSMs with their own converter
        if 'search_hit' not in psm.keys():
            nohit_PSM_to_SAM(psm,file)
        else:
            for row in psm['search_hit']:
                if rm_duplicates=="Y":
                    dup={}
                for p in range(0,len(row['proteins'])):
                    decoy=0
                    # convert decoys with decoy-specific convertor
                    for d in decoy_annotation:
                        if d in row['proteins'][p]['protein'].upper():
                            decoy=1
                            key=row['proteins'][p]['protein'].upper().split(d)[1]
                            if key in transcript_hash.keys():
                                temp_result= decoy_PSM_to_SAM(psm,row,key,transcript_hash,exon_hash,allowed_mismatches,
                                                 map_decoy)
                                if rm_duplicates=="Y":
                                    dup_key= str(temp_result[9])+"_"+str(temp_result[2])+"_"+"_"+str(temp_result[3])
                                    if dup_key not in dup.keys():
                                        dup[dup_key]=1
                                        write_psm(temp_result,file)
                                else:
                                    write_psm(temp_result,file)
                            else:
                                temp_result=unannotated_PSM_to_SAM(psm,row,decoy)
                                if rm_duplicates=="Y":
                                    dup_key= str(temp_result[9])+"_"+str(temp_result[2])+"_"+"_"+str(temp_result[3])
                                    if dup_key not in dup.keys():
                                        dup[dup_key]=1
                                        write_psm(temp_result,file)
                                else:
                                    write_psm(temp_result,file)

                    if decoy==0:
                        key=row['proteins'][p]['protein']
                        # Filter out PSM where transcript sequences were not found/ non-existent
                        if key not in transcript_hash.keys():
                            write_psm(unannotated_PSM_to_SAM(psm,row,decoy),file)
                        else:
                            if three_frame_translation=="Y":
                                protein_hit=map_peptide_to_protein_3frame(row['peptide'],transcript_hash[key]['transcript_seq'],allowed_mismatches,transcript_hash[key]['strand'])
                            else:
                                protein_hit=map_peptide_to_protein(row['peptide'],transcript_hash[key]['protein_seq']
                                                                   ,allowed_mismatches)
                            if len(protein_hit)==0:
                                write_psm(unannotated_PSM_to_SAM(psm,row,decoy),file)
                            else:
                                # map peptide on protein and retrieve hit position, iterate over all hits
                                for phit in protein_hit:
                                    temp_result=[None]*23
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
                                                                 transcript_hash[key]['strand'],row['peptide'])
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
                                        temp_result[9]=_reverse_complement_(str(transcript_hash[key]['transcript_seq']\
                                                   [phit_loc:(phit_loc+(len(row['peptide'])*3))]))
                                    #QUAL
                                    temp_result[10]='*'
                                    #
                                    #Mandatory proteomics specific columns added to the proBam format
                                    #
                                    #NH: number of genomic location the peptide mapping to
                                    temp_result[11]='NH:i:'+str(len(row['proteins'])+len(phit)-1)
                                    #XL: number of peptides the spectrum mapping to
                                    temp_result[12]='XL:i:'+str(len(psm['search_hit']))
                                    #XP; peptide sequence
                                    temp_result[13]='XP:Z:'+row['modified_peptide']
                                    #XR: reference peptide sequence
                                    temp_result[14]='XR:Z:'+row['peptide']
                                    #XS: PSM score
                                    temp_result[15]=create_XS(row['search_score'])
                                    #XQ: PSM-Qvalue
                                    temp_result[16]='XQ:f:'+'NA'
                                    #XC: Peptide charge
                                    temp_result[17]='XC:i:'+str(psm['assumed_charge'])
                                    #XA: Whether the peptide is well annotated
                                    temp_result[18]=create_XA(phit[1])
                                    #XM: Modification
                                    temp_result[19]=create_XM(row['modifications'])
                                    #XN: number of mis-cleavages
                                    if 'num_missed_cleaveges' in row.keys():
                                        temp_result[20]='XN:i:'+str(row['num_missed_cleavages'])
                                    else:
                                        temp_result[20]='XN:i:0'
                                    #XT: non/semi/tryptic
                                    temp_result[21]="XT:i:NA"
                                    #XG: Petide type
                                    temp_result[22]=create_XG(phit[1])
                                    # remove duplicates if rm_duplicates=Y
                                    if rm_duplicates=="Y":
                                        dup_key= dup_key= str(temp_result[9])+"_"+\
                                                          str(temp_result[2])+"_"+str(temp_result[3])
                                        if dup_key not in dup.keys():
                                            dup[dup_key]=1
                                            write_psm(temp_result,file)
                                    else:
                                        write_psm(temp_result,file)
    file.close()


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
# Function to determine peptide annotation and create XA string
#
def create_XA(mismatch):
    '''
    :param mismatch: number of allowed mismatches
    :return: XA
    '''
    XA='XA:i:'
    if mismatch==0:
        XA+='0'
    else:
        XA+='1'
    return XA
#
# Function to determine PSM-score and create XS string
#
def create_XS(score):
    '''
    :param score: psm score
    :return: xcorr score ( if existing)
    '''
    XS='XS:f:'
    if 'xcorr' in score.keys():
        XS+=str(score['xcorr'])
    else:
        XS+='NA'
    return XS
#
# Function to analyze peptide type and create XG string
#
def create_XG(hamming_distance):
    '''
    :param hamming_distance: hamming distance to use
    :return: the hamming distance
    '''
    XG='XG:Z:'
    if hamming_distance==0:
        XG+='N'
    else:
        XG+='V'
    return XG
#
# Function to analyze modification and create XA-string
#
def create_XM(modifications):
    '''
    :param modifications: psm modification
    :return: modification in proBAM format
    '''
    XM='XM:Z:'
    if modifications==[]:
        XM+='NA'
    else:
        for mod in modifications:
            if ';' in XM:
                XM+=';'
            XM+=str(mod['position'])
            XM+=';'
            XM+=str(mod['mass'])
    return XM
#
# Function to calculate FLAG of PSM
#
def calculate_FLAG(strand,rank,decoy):
    '''
    :param strand: transcript strand
    :param rank: psm rank
    :param decoy: boolean decoy
    :return: FLAG
    '''
    strand=int(strand)
    FLAG=0
    if strand==-1:
        FLAG=FLAG+16
    if rank!=1:
        FLAG=FLAG+256
    if decoy==1:
        FLAG=FLAG+1024
    return FLAG

#
# Function that maps peptide on the corresponding protein, alowing mismatches ( as specified)
# and returns the first left base position mapped
#
def map_peptide_to_protein(peptide_seq,protein_seq,allowed_mismatches):
    '''
    :param peptide_seq: peptide sequence
    :param protein_seq: protein sequence
    :param allowed_mismatches: allowed mismatches
    :return: location of mapped peptide
    '''
    hits=[]
    pep_length=len(peptide_seq)
    prot_length=len(protein_seq)
    for i in range(0,prot_length-pep_length):
        if hamming(peptide_seq,protein_seq[i:pep_length+i]) <= allowed_mismatches:
            hits.append([i,hamming(peptide_seq,protein_seq[i:pep_length+i])])
    return hits

#
#Function that maps peptide on the corresponding protein, alowing mismatches ( as specified)
# and returns the first left base position mapped after translating the transcript sequence in 3-frames
#
def map_peptide_to_protein_3frame(peptide_seq,transcript_seq,allowed_mismatches,strand):
    size_adjust=-1    # adjust size of transcript for starting at +1/+2frame
    hits=[]
    pep_length=len(peptide_seq)
    frame=[0]*3
    frame[0]=standard_code.translate(transcript_seq)
    frame[1]=standard_code.translate(transcript_seq[1:])
    frame[2]=standard_code.translate(transcript_seq[2:])
    for f in frame:
        size_adjust+=1
        for i in range(0,(len(f)-pep_length)):
            if hamming(peptide_seq,f[i:pep_length+i]) <= allowed_mismatches:
                adjusted_hit_pos=(i*3)+size_adjust
                hits.append([adjusted_hit_pos,hamming(peptide_seq,f[i:pep_length+i])])
    return hits

#
# Function returns genomic position of the leftmost base
#
#
def calculate_genome_position(phit,strand,offset,start_exon_rank,peptide,exons,chr,three_frame_translation):
    '''
    :param phit: location of peptide hit
    :param strand: transcript strand
    :param offset: CDS offset
    :param start_exon_rank: rank of first exon containing the CDS
    :param peptide: peptide sequence
    :param exons: exon hash
    :param chr: transcript chr
    :return: genomic start position
    '''
    if three_frame_translation=='Y':
        if strand=='1':
            tr_pos=(phit)
            pointer=0
            #start_exon_rank=0
        else:
            tr_pos=(phit)+(len(peptide)*3)-1
            #start_exon_rank=len(exons)+1
            pointer=0
        gen_pos=0
        offset=0
        start_exon_rank=1

    else:
        if strand=='1':
            tr_pos=(phit*3)
        else:
            tr_pos=(phit*3)+(len(peptide)*3)-1
        pointer=1
        gen_pos=0

    # iterate over exons till peptide position is encoutered
    exons=sorted(exons,key=lambda x: int(x[2]))
    for exon in exons:
        exon=[int(numeric_string) for numeric_string in exon]
        if exon[2]<start_exon_rank:
            continue
        elif exon[2]==start_exon_rank:
            if strand=='1':
                exon[0]=exon[0]+offset
            else:
                exon[1]=exon[1]-offset
        if tr_pos>((exon[1]-exon[0])+pointer):
            pointer=pointer+(exon[1]-exon[0]+1)
            continue
        else:
            remember_exon=exon
            if strand=='1':
                gen_pos=exon[0]+((tr_pos-pointer))
            else:
                gen_pos=exon[1]-((tr_pos-pointer))
            break

    # below print tests are generated to get the left most basepare
    '''
    print '*************************************************'
    print 'test:',start_exon_rank,phit,strand
    print pointer,tr_pos, offset
    print remember_exon
    print gen_pos,peptide,strand,chr
    print '*************************************************'
    print ''
    '''
    return str(gen_pos)

#
# Function to compute the CIGAR string of a protein
#
def compute_cigar(gen_pos,exons,strand,peptide):
    '''
    :param gen_pos: genomic start position
    :param exons: exon dictionairy
    :param strand: transcript strand
    :param peptide: peptide sequence
    :return: CIGAR string
    '''
    gen_pos=int(gen_pos)
    hit=0
    cigar=''
    pointer=0
    debug=0
    # debug section for wrong CIGAR
    #if peptide=="TLISDIEAVK":
    #    debug=1
    length= int(len(peptide)*3)
    exons=sorted(exons,key=lambda x: int(x[2]))
    if strand=='-1':
        adjusted_exons=list(reversed(exons))
    else:
        adjusted_exons=exons
    for exon in adjusted_exons:
        exon=[int(numeric_string) for numeric_string in exon]
        if hit ==0:
            if gen_pos >= exon[1]:
                continue
            else:
                if (gen_pos+length)<=exon[1]+1:
                    cigar=str(length)+'M'
                    break
                else:
                    pointer=(exon[1]-gen_pos+1)
                    cigar=cigar+(str(pointer)+'M')
                    length=length-pointer
                    prev_exon_end=exon[1]
                    prev_rank=exon[2]
                    hit=1
                    continue
        if hit!=0:
            cigar+=(str((exon[0]-prev_exon_end-1))+'N')
            if (exon[0]+(length)-1)>exon[1]+1:
                cigar=cigar+(str((exon[1]-exon[0]))+'M')
                length=length-(exon[1]-exon[0])
                continue
            else:
                cigar=cigar+(str(length)+'M')
                break
    #debug section fro wrong CIGAR
    '''if debug==1:
        print strand
        print peptide
        print cigar
        print exons
        print gen_pos'''
    return cigar

#
# Function to convert unmapped PSM to SAM
#

def nohit_PSM_to_SAM(psm,file):
    '''
    :param psm: psm file
    :param file: output file
    :return: converted psm to SAM if no hits found
    '''
    temp_result=[None]*23
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
    #XL: number of peptides the spectrum mapping to
    temp_result[12]='XL:i:0'
    #XP; peptide sequence
    temp_result[13]='XP:Z:NA'
    #XR: reference peptide sequence
    temp_result[14]='XR:Z:NA'
    #XS: PSM score
    temp_result[15]='XS:f:0'
    #XQ: PSM-Qvalue
    temp_result[16]='XQ:f:0'
    #XC: Peptide charge
    temp_result[17]='XC:i:0'
    #XA: Whether the peptide is well annotated
    temp_result[18]='XA:i:0'
    #XM: Modification
    temp_result[19]='XM:Z:NA'
    #XN: number of mis-cleavages
    temp_result[20]='XN:i:0'
    #XT: non/semi/tryptic
    temp_result[21]="XT:i:0"
    #XG: Petide type
    temp_result[22]="XG:Z:U"
    write_psm(temp_result,file)


#
# Function to convert unannotated PSMs to SAM
#
def unannotated_PSM_to_SAM(psm,row,decoy):
    '''
    :param psm: psm dictionairy
    :param row: unnanotated PSM row
    :param decoy: decoy boolean
    :param file: output file
    :return: sam of unnanotated PSM
    '''
    decoy=int(decoy)
    temp_result=[None]*23
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
    #XL: number of peptides the spectrum mapping to
    temp_result[12]='XL:i:'+str(len(psm['search_hit']))
    #XP; peptide sequence
    temp_result[13]='XP:Z:'+row['modified_peptide']
    #XR: reference peptide sequence
    temp_result[14]='XR:Z:'+row['peptide']
    #XS: PSM score
    temp_result[15]=create_XS(row['search_score'])
    #XQ: PSM-Qvalue
    temp_result[16]='XQ:f:0'
    #XC: Peptide charge
    temp_result[17]='XC:i:'+str(psm['assumed_charge'])
    #XA: Whether the peptide is well annotated
    temp_result[18]='XA:i:2'
    #XM: Modification
    temp_result[19]=create_XM(row['modifications'])
    #XN: number of mis-cleavages
    if 'num_missed_cleaveges' in row.keys():
        temp_result[20]='XN:i:'+str(row['num_missed_cleavages'])
    else:
        temp_result[20]='XN:i:0'
    #XT: non/semi/tryptic
    temp_result[21]="XT:i:0"
    #XG: Petide type
    if decoy==1:
        temp_result[22]="XG:Z:U"
    else:
        temp_result[22]="XG:Z:D"
    return temp_result

#
# Function to convert decoy PSM to SAM format
#

def decoy_PSM_to_SAM(psm,row,key,transcript_hash,exon_hash,allowed_mismatches,map_decoy):
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
        protein_hit=map_peptide_to_protein(row['peptide'][::-1],transcript_hash[key]['protein_seq'],allowed_mismatches)
    else:
        protein_hit=[]
    if len(protein_hit)==0:
         return unannotated_PSM_to_SAM(psm,row,1)
    else:
        # map peptide on protein and retrieve hit position, iterate over all hits
        for phit in protein_hit:
            temp_result=[None]*23
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
                temp_result[9]=_reverse_complement_(str(transcript_hash[key]['transcript_seq']\
                           [(phit[0]*3):((phit[0]*3)+(len(row['peptide'])*3))]))
            #QUAL
            temp_result[10]='*'
            #
            #Mandatory proteomics specific columns added to the proBam format
            #
            #NH: number of genomic location the peptide mapping to
            temp_result[11]='NH:i:'+str(len(row['proteins'])+len(phit)-1)
            #XL: number of peptides the spectrum mapping to
            temp_result[12]='XL:i:'+str(len(psm['search_hit']))
            #XP; peptide sequence
            temp_result[13]='XP:Z:'+row['modified_peptide']
            #XR: reference peptide sequence
            temp_result[14]='XR:Z:'+row['peptide']
            #XS: PSM score
            temp_result[15]=create_XS(row['search_score'])
            #XQ: PSM-Qvalue
            temp_result[16]='XQ:f:'+'NA'
            #XC: Peptide charge
            temp_result[17]='XC:i:'+str(psm['assumed_charge'])
            #XA: Whether the peptide is well annotated
            temp_result[18]=create_XA(phit[1])
            #XM: Modification
            temp_result[19]=create_XM(row['modifications'])
            #XN: number of mis-cleavages
            if 'num_missed_cleaveges' in row.keys():
                temp_result[20]='XN:i:'+str(row['num_missed_cleavages'])
            else:
                temp_result[20]='XN:i:0'
            #XT: non/semi/tryptic
            temp_result[21]="XT:i:NA"
            #XG: Petide type
            temp_result[22]='XG:Z:D'
            return temp_result
#
# Function to calculate distance between 2 strings
#
def hamming(str1, str2):
    '''
    :param str1: first string
    :param str2: second string
    :return: the hamming distance
    '''
    assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))

#
# Create SAM header
#
def create_SAM_header(file,version,database,sorting_order,database_v,species):
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
    header.append('@PG\tID:proBamPy\tVN:1.0')
    header.append('@GA\tAS:'+str(database)+'\tVN:'+str(database_v))
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
    pysam.sort((directory+name+'.bam'),(directory+name+'.sorted'))
    pysam.index(directory+name+'.sorted.bam')

#
# reverse complements a DNA string
#
def _reverse_complement_(DNA):
    temp=DNA[::-1]
    DNA=""
    for i in temp:
        if i=="A":
            DNA+='T'
        elif i=="T":
            DNA+='A'
        elif i=="G":
            DNA+='C'
        else:
            DNA+='G'
    return DNA

####################
### MAIN PROGRAM ###
####################

if __name__=='__main__':
    get_input_variables()
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
    create_SAM_header(file,version,database,sorting_order,database_v,species)
    PSM2SAM(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,map_decoy,rm_duplicates,
            three_frame_translation)
    sam_2_bam(directory,name)



    print("proBAM conversion succesful")
    print("%f seconds" % (time.time() - start_time))         # output script run time
