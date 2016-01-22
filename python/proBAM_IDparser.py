__author__="vladie"
import proBAM_ENSEMBL
import proBAM_biomart
import re
import sys

#
# Find accession encryption pattern and database
#

def parseID(psm_hash,species,database,decoy_annotation,database_v):
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
    if len(PID_list)<100:
        max=len(PID_list)
    else:
        max=100
    # parse | and _ out of identifier and put in a list of terms
    while found==0 and count<max:
        if "|" in PID_list[count]:
            decrypted_list=PID_list[count].split("|")
            find_result=find_annotation(decrypted_list,species)
            if find_result[0]==1:
                found=1
                break
        elif "_" in PID_list[count]:
            decrypted_list=PID_list[count].split("_")
            find_result=find_annotation(decrypted_list,species)
            if find_result[0]==1:
                found=1
                break
        else:
            decrypted_list=[PID_list[count]]
            find_result=find_annotation(decrypted_list,species)
            if find_result[0]==1:
                found=1
                break
        count+=1
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
                        psm['proteins'][i]['protein']=update_protein_accession(psm['proteins'][i]['protein'],
                                                                               decoy_annotation,find_result[2])
                        hit=0
                        for d in decoy_annotation:
                            if d in psm['proteins'][i]['protein'].upper():
                                ID=psm['proteins'][i]['protein'].upper().split(d)[1]
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
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(protein_ID,'transcript',database_v,species)
            elif find_result[1]=="ENSEMBL_PR":
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(protein_ID,'protein',database_v,species)
            elif find_result[1]=="UNIPROT":
                conversion=id_map('UNIPROT/SWISSPROT','ENSEMBL',protein_ID,psm_hash,species,decoy_annotation,database_v)
                psm_protein_id=conversion[0]
                psm_hash=conversion[1]
                annotation= proBAM_ENSEMBL.prepareAnnotationENSEMBL(psm_protein_id,'transcript',database_v,species)
            return [psm_hash,annotation]
        # raise error if database unsupported
        else:
            raise ValueError('Currently supported annotation databases: \n ENSEMBL')



#
# Find annotation in list of terms
#

def find_annotation(term_list,species):
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
            if ensembl_prefix[0][0] in term_list[pos]:
                print "Identified ENSEMBL transcript IDs"
                result=[1,"ENSEMBL_TR",pos]
                found=1
                break
            elif ensembl_prefix[0][1] in term_list[pos]:
                print "Identified ENSEMBL protein IDs"
                result=[1,"ENSEMBL_PR",pos]
                found=1
                break
            elif re.match("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",term_list[pos]) \
                            is not None:
                print " Identified UNIPROT/SWISSPROT IDs"
                result=[1,"UNIPROT",pos]
                found=1
                break

    return result



#
# Parse protein accesions
#

def update_protein_accession(accession,decoy_annotation,hit):
    decoy=0
    if any(decoy in accession.upper() for decoy in decoy_annotation):
        decoy=1
    if '|' in accession:
        accession=accession.split('|')[hit]
    if "_" in accession:
        if decoy==1:
            if "_REVERSED" in accession:
                accession="DECOY_"+accession.split('_')[hit]
            else:
                accession="DECOY_"+accession.split('_')[hit+1]
        else:
            accession=accession.split('_')[hit]
    return accession

#
# Retrieve matched protein IDs
#
def get_protein_ID(psm_hash,decoy_annotation):
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
                            ID=psm['proteins'][i]['protein'].upper().split(d)[1]
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
def id_map(from_annotation,to_annotation,psm_protein_id,psm_hash,species,decoy_annotation,database_v):
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
            mapped_id=proBAM_biomart.id_map_ensembl("uniprot_swissprot",database_v,species,psm_protein_id)
            for row in mapped_id:
                row=row.split("\t")
                if row[0]!="":
                    if row[2] in map.keys():
                        if map[row[2]][1]<row[1]:
                            map[row[2]]=row
                    else:
                        map[row[2]]=row
    for key in psm_hash:
        if 'search_hit' in key.keys():
            for psm in key['search_hit']:
                for i in range (0,len(psm['proteins'])):
                    hit=0
                    for d in decoy_annotation:
                        if d in psm['proteins'][i]['protein'].upper():
                            ID=psm['proteins'][i]['protein'].upper().split(d)[1]
                            if ID in map.keys():
                                psm['proteins'][i]['protein']="DECOY_"+map[ID][0]
                                new_psm_protein_id.append(map[ID][0])
                            hit=1
                    if hit==0:
                        if psm['proteins'][i]['protein'] in map.keys():
                            new_psm_protein_id.append(map[psm['proteins'][i]['protein']][0])
                            psm['proteins'][i]['protein']=map[psm['proteins'][i]['protein']][0]


    #retain only unique IDs
    return [new_psm_protein_id,psm_hash]

#
# When an ID maps to multiple transcript IDs, fetch the longest transcript
#
#TODO after conformation, retrieve only the longest sequences (bp or aa ? )
def fetch_longest_transcript(ID):
    return ID[0]
