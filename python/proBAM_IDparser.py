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

__author__="vladie"
import proBAM_ENSEMBL
import proBAM_biomart
import re
import sys
from fnmatch import fnmatch

#
# Find accession encryption pattern and database
#

def parseID(psm_hash,species,database,decoy_annotation,database_v,three_frame_translation,pre_picked_annotation):
    '''
    :param psm_hash: psm hash (see: proBAM_input.get_PSM_hash())
    :param species: species
    :param database: annotation database name
    :param decoy_annotation: decoy annotation list
    :param database_v: annotation database version
    :return: [PSM_hash_with_updated_identifiers,[transcript_annotation_hash,exon_hash]]
    '''
    # create list with IDs, for faster access
    PID_list=[]
    for key in psm_hash:
            if 'search_hit' in key.keys():
                for psm in key['search_hit']:
                    for i in range (0,len(psm['proteins'])):
                       PID_list.append(psm['proteins'][i]['protein'])
    found=0
    count=0
    # set maximum amount of identifiers to be checked before giving up identifier identification
    if len(PID_list)<1000:
        max=len(PID_list)
    else:
        max=1000
    # parse | and _ out of identifier and put in a list of terms
    while found==0 and count<max:
        if "|" in PID_list[count]:
            decrypted_list=PID_list[count].split("|")
            find_result=_find_annotation_(decrypted_list,species,pre_picked_annotation)
            if find_result[0]==1:
                found=1
                break
        elif "_" in PID_list[count]:
            decrypted_list=PID_list[count].split("_")
            find_result=_find_annotation_(decrypted_list,species,pre_picked_annotation)
            if find_result[0]==1:
                found=1
                break
        else:
            decrypted_list=[PID_list[count]]
            find_result=_find_annotation_(decrypted_list,species,pre_picked_annotation)
            if find_result[0]==1:
                found=1
                break
        count+=1
    # print PID_list[0:100]
    # raise error max amount of identifiers checked and nog identifiers identified
    if found==0:
        raise ValueError('Protein ID annotations not recognized\n'\
            'currently supported Protein ID annotations:\n'\
            'ENSEMBL protein identifiers\n ENSEMBL transcript identifiers\n UNIPROT/SWISSPROT identifiers ')
        sys.exit()

    # if identifiers identified, update PSM_hash to the right identifiers
    elif found==1:
        protein_ID=[]
        for key in psm_hash:
            if 'search_hit' in key.keys():
                for psm in key['search_hit']:
                    for i in range (0,len(psm['proteins'])):
                        psm['proteins'][i]['protein']=_update_protein_accession_(psm['proteins'][i]['protein'],
                                                                               decoy_annotation,find_result[2])
                        hit=0
                        for d in decoy_annotation:
                            if d in psm['proteins'][i]['protein'].upper():
                                ID=psm['proteins'][i]['protein'].upper().strip(d)
                                protein_ID.append(ID)
                                hit=1
                        if hit==0:
                            protein_ID.append(psm['proteins'][i]['protein'])


        #retain only unique IDs
        unique=set(protein_ID)
        protein_ID=list(unique)
        # if database == ENSEMBL use correct function to retrieve anntotation information
        if database=="ENSEMBL":
            if find_result[1]=="ENSEMBL_TR":
                conversion = _id_map_('ENSEMBL', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,database_v)
                id_map = conversion[0]
                psm_hash = conversion[1]
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(id_map.values(),'transcript',database_v,species,three_frame_translation)
            elif find_result[1]=="ENSEMBL_PR":
                conversion = _id_map_('ENSEMBL', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,database_v)
                id_map = conversion[0]
                psm_hash = conversion[1]
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(id_map.values(),'protein',database_v,species,three_frame_translation)
            elif find_result[1]=="UNIPROT":
                conversion=_id_map_('UNIPROT/SWISSPROT','ENSEMBL',protein_ID,psm_hash,species,decoy_annotation,database_v)
                id_map=conversion[0]
                psm_hash=conversion[1]
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(id_map.values(),'transcript',database_v,species,three_frame_translation)
            return [psm_hash,annotation,id_map]
        # raise error if database unsupported
        else:
            raise ValueError('Currently supported annotation databases: \n ENSEMBL')



#
# Find annotation in list of terms
#

def _find_annotation_(term_list,species,pre_picked_annotation):
    '''
    :param term_list: list of terms of possible identifiers
    :param species: species
    :return: list [found_boolean,Type_identified_identifier,position of identifier term]
    '''
    ensembl_prefix=proBAM_ENSEMBL.get_Ensembl_prefix(species)
    found=0
    result=[0,0,0]
    for pos in range(0,len(term_list)):
        if found==0:
            if ensembl_prefix[0][0] in term_list[pos] and \
                    (pre_picked_annotation=="First" or pre_picked_annotation=="Ensembl_tr"):
                print "Identified ENSEMBL transcript IDs"
                result=[1,"ENSEMBL_TR",pos]
                found=1
                break
            elif ensembl_prefix[0][1] in term_list[pos] and \
                (pre_picked_annotation == "First" or pre_picked_annotation == "Ensembl_pr"):
                print "Identified ENSEMBL protein IDs"
                result=[1,"ENSEMBL_PR",pos]
                found=1
                break
            elif (pre_picked_annotation=="First" or pre_picked_annotation=="UniProt"):
                if re.match(".*[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}.*",term_list[pos]) \
                            is not None:
                    print " Identified UNIPROT/SWISSPROT IDs"
                    result=[1,"UNIPROT",pos]
                    found=1
                    break
    return result



#
# Parse protein accesions
#

def _update_protein_accession_(accession,decoy_annotation,hit):
    '''
    :param accession: list of possible protein ID terms
    :param decoy_annotation: list of decoy annotations
    :param hit: position of protein ID
    :return: protein ID
    '''
    is_decoy=0
    if any(decoy in accession.upper() for decoy in decoy_annotation):
        is_decoy=1
    forward_decoys=[]
    backward_decoys=[]
    for decoy in decoy_annotation:
        if fnmatch(decoy,"_*"):
            backward_decoys.append(decoy.replace('_',''))
        elif fnmatch(decoy,"*_"):
            forward_decoys.append(decoy.replace('_',''))
    if '|' in accession:
        accession=accession.split('|')[hit]
    if "_" in accession:
        if is_decoy==1:
            for i in range(0,len(accession.split('_'))):
                if any(b_decoy in accession.split('_')[i] for b_decoy in backward_decoys):
                    accession="DECOY_"+accession.split('_')[i-1]
                elif any(f_decoy in accession[i].split('_') for f_decoy in forward_decoys):
                    accession="DECOY_"+accession.split('_')[i+1]
        else:
            accession=accession.split('_')[hit]
    return accession

#
# Retrieve matched protein IDs
#
def _get_protein_ID_(psm_hash,decoy_annotation):
    '''
    :param psm_hash: dictionairy containing PSM file
    :param decoy_annotation: list of decoy annotation tags
    :return: list of unique protein IDs
    '''

    #retrieve protein IDs from PSM HASH
    protein_ID=[]
    for key in psm_hash:
        if 'search_hit' in key.keys():
            for psm in key['search_hit']:
                for i in range (0,len(psm['proteins'])):
                    hit=0
                    for d in decoy_annotation:
                        if d in psm['proteins'][i]['protein'].upper():
                            ID=psm['proteins'][i]['protein'].upper().strip(d)
                            protein_ID.append(ID)
                            hit=1
                    if hit==0:
                        protein_ID.append(psm['proteins'][i]['protein'])
    #retain only unique IDs
    unique=set(protein_ID)
    #print unique
    return list(unique)

#
# Function to map protein IDs between annotations
#
def _id_map_(from_annotation,to_annotation,psm_protein_id,psm_hash,species,decoy_annotation,database_v):
    '''
    :param from_annotation: supplied annotation (i.e. swissprot)
    :param to_annotation: target annotation (i.e. ENSEMBL)
    :param psm_protein_id: list of protein IDS
    :param psm_hash: dictionairy of protein IDs mapped onto ENSEMBL
    :param species: species name
    :param decoy_annotation: list of decoy annotations
    :param database_v: database version
    :return: dictionairy of protein ID coversion
    '''
    #psm_hash.reset()
    new_psm_protein_id=[]

    print "Commencing ID conversion from "+str(from_annotation)+" to "+str(to_annotation)
    map={}
    if to_annotation=="ENSEMBL":
        if from_annotation=="UNIPROT/SWISSPROT":
            temp_map={}
            mapped_id=proBAM_biomart.id_map_ensembl("uniprot_swissprot",database_v,species,psm_protein_id)
            for row in mapped_id:
                row=row.split("\t")
                if row[0]!="":
                    if row[2] in temp_map.keys():
                        if temp_map[row[2]][1]<row[1]:
                            temp_map[row[2]]=row
                    else:
                        temp_map[row[2]]=row
            for key in temp_map:
                map[key]=temp_map[key][0]

        if from_annotation=="ENSEMBL":
            for id in psm_protein_id:
                map[id]=id.split('.')[0]


    return [map,psm_hash]

#
# When an ID maps to multiple transcript IDs, fetch the longest transcript
#
#TODO after conformation, retrieve only the longest sequences (bp or aa ? )
def fetch_longest_transcript(ID):
    return ID[0]
