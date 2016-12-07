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
import urllib,urllib2
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
    decryption=[]
    while found==0 and count<max:
        find_result = _find_annotation_(PID_list[count],species, pre_picked_annotation)
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
                                continue
                                hit=1
                        if hit==0:
                            protein_ID.append(psm['proteins'][i]['protein'])

        # if database == ENSEMBL use correct function to retrieve anntotation information
        if database=="ENSEMBL":
            if find_result[1]=="ENSEMBL_TR":
                id_map = _id_map_('ENSEMBL', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,database_v)
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),'transcript',
                                                                    database_v,species,three_frame_translation)
            elif find_result[1]=="ENSEMBL_PR":
                id_map = _id_map_('ENSEMBL', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,database_v)
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),'protein',
                                                                    database_v,species,three_frame_translation)
            elif find_result[1]=="UNIPROT":
                id_map = _id_map_('UNIPROT','ENSEMBL',protein_ID,psm_hash,species,decoy_annotation,database_v)
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map),'transcript',
                                                                    database_v,species,three_frame_translation)
            elif find_result[1] == "UNIPROT_ENTRY":
                id_map = _id_map_('UNIPROT_ENTRY', 'ENSEMBL', protein_ID, psm_hash, species, decoy_annotation,
                                      database_v)
                annotation = proBAM_ENSEMBL.prepareAnnotationENSEMBL(_get_transcript_ids_from_map(id_map), 'transcript',
                                                                     database_v, species, three_frame_translation, )
            return [psm_hash,annotation,id_map]
        # raise error if database unsupported
        else:
            raise ValueError('Currently supported annotation databases: \n ENSEMBL')

#
# Create a list of unique transcript ID's from MAP
#
def _get_transcript_ids_from_map(map):
    transcript_ids=[]
    for value in map.values():
        transcript_ids+=value
    transcript_ids=list(set(transcript_ids))
    return transcript_ids

#
# Find annotation in list of terms
#

def _find_annotation_(accession,species,pre_picked_annotation):
    '''
    :param term_list: list of terms of possible identifiers
    :param species: species
    :return: list [found_boolean,Type_identified_identifier,position of identifier term]
    '''
    found=0
    ensembl_prefix=proBAM_ENSEMBL.get_Ensembl_prefix(species)
    result=[0,'']
    if re.findall(ensembl_prefix[0][0]+"[0-9]*",accession)!=[] and \
            (pre_picked_annotation=="First" or pre_picked_annotation=="Ensembl_tr") and found==0:
        print "Identified ENSEMBL transcript IDs"
        result=[1,"ENSEMBL_TR",ensembl_prefix[0][0]+"[0-9]*"]
        found=1
    if re.findall(ensembl_prefix[0][1]+"[0-9]*",accession)!=[] and \
        (pre_picked_annotation == "First" or pre_picked_annotation == "Ensembl_pr") and found==0 :
        print "Identified ENSEMBL protein IDs"
        result=[1,"ENSEMBL_PR",ensembl_prefix[0][1]+"[0-9]*"]
        found=1
    if (pre_picked_annotation=="First" or pre_picked_annotation=="UniProt_ACC") and found==0:
        if re.findall("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",accession)!=[]:
            print " Identified UNIPROT IDs"
            result=[1,"UNIPROT","[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]"]
            found=1
    if (pre_picked_annotation == "First" or pre_picked_annotation == "UniProt_Entry") and found == 0:
        if re.findall("[A-Z0-9]{1,10}"+"_"+_get_uniprot_postfix_(species),accession)!=[]:
            print "Identified UNIPROT Entry Names"
            result=[1,"UNIPROT_ENTRY","[A-Z0-9]{1,10}"+"_"+_get_uniprot_postfix_(species)]
            found=1
    return result



#
# Parse protein accesions
#

def _update_protein_accession_(accession,decoy_annotation,regex):
    '''
    :param accession: list of possible protein ID terms
    :param decoy_annotation: list of decoy annotations
    :param hit: position of protein ID
    :return: protein ID
    '''
    is_decoy=0

    if any(decoy in accession.upper().replace('-','_') for decoy in decoy_annotation):
        is_decoy=1
    new_accession=re.findall(regex,accession)
    if new_accession!=[] and new_accession[0]!="P00000":
        accession= new_accession[0]
        if is_decoy==1:
            accession="DECOY_"+accession
    return accession
#
# Get UNIPROT entry names postfixes
#
def _get_uniprot_postfix_(species):
    if species=='homo_sapiens':
        return "HUMAN"
    elif species=='mus_musculus':
        return "MOUSE"
    elif species=='drosophila_melanogaster':
        return "DROME"
    elif species=='danio_rerio':
        return "DANRE"
    elif species=='arabidopsis_thaliana':
        return "ARATH"
    else:
        raise ValueError('Species not recognized')
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
        if from_annotation=="UNIPROT":
            temp_map={}
            mapped_id=proBAM_biomart.id_map_ensembl("uniprot_swissprot",database_v,species,psm_protein_id)
            for row in mapped_id:
                if row[0]!="":
                    if row[2] in temp_map.keys():
                            temp_map[row[2]].append(row)
                    else:
                        temp_map[row[2]]=[row]
            for key in temp_map:
                map[key]=temp_map[key][0]

        if from_annotation=="UNIPROT_ENTRY":
            from bioservices import UniProt
            u=UniProt()
            to_translate=[]
            for id in psm_protein_id:
                if re.findall("[A-Z0-9]{1,10}" + "_" + _get_uniprot_postfix_(species), id) != []:
                    to_translate.append(re.findall("[A-Z0-9]{1,10}" + "_" + _get_uniprot_postfix_(species), id)[0])

            # map uniprot_entries to up-to-date accession
            to_translate=list(set(to_translate))
            accession_update_hash=u.mapping('ACC+ID','ACC',to_translate)
            for entry in to_translate:
                if entry not in accession_update_hash:
                    try:
                        accession_update_hash[entry]=[get_updated_entry_name(entry)]
                    except urllib2.HTTPError:
                        pass
            to_translate=[]

            #map accession to Ensembl transcript
            for value in accession_update_hash.values():
                to_translate+=value
            to_translate=list(set(to_translate))

            # remap Ensembl transcript to uniprot entries
            temp_map=u.mapping('ACC','ENSEMBL_TRS_ID',to_translate)
            for accession in accession_update_hash:
                for i in accession_update_hash[accession]:
                    if i in temp_map:
                        if accession in map:
                            map[accession]+=temp_map[i]
                        else:
                            map[accession]=temp_map[i]

        if from_annotation=="ENSEMBL":
            for id in psm_protein_id:
                map[id]=[id]
    return map
#
# given a UniProt entry name, retrieves the stable accession ID
#
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

