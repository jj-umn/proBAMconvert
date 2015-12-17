__author__ = 'vladie'

from pyteomics import mzid,mass
import re

#
# Parse mzid file into dictionairy
#

def get_PSM_mzid(psm_file):
    '''
    :param psm_file: mzid file
    :return: dictionairy of parsed mzid file, suitable for proBAMconvert
    '''
    PSM=mzid.read(psm_file)
    psm_hash=[]
    for row in PSM:
        #print row
        temp_hash={"assumed_charge":row['SpectrumIdentificationItem'][0]['chargeState'],"spectrum":row['id'],"search_hit":[]}
        for psm in row["SpectrumIdentificationItem"]:
            proteins=[]
            for protein in psm["PeptideEvidenceRef"]:
                if '|' in protein['peptideEvidence_ref']:
                    if len(protein['peptideEvidence_ref'].split("|")[1])!=1:
                        proteins.append({"protein":protein['peptideEvidence_ref'].split("|")[1]})
                    else:
                        proteins.append({"protein":protein['peptideEvidence_ref'].split("|")[0]})
                if '_' in protein['peptideEvidence_ref']:
                    proteins.append({"protein":protein['peptideEvidence_ref'].split("_")[1]})
            peptide=_get_peptide_sequence_(psm['peptide_ref'])
            modifications=_get_peptide_modifications_(psm['peptide_ref'])
            score=_get_score_(psm)

            temp_hash['search_hit'].append({"hit_rank":psm['rank'],"modifications":modifications,
                                            "modified_peptide":psm['peptide_ref'],"peptide":peptide,
                                            "search_score":score,"proteins":proteins,"num_missed_cleavages":"0"})
        #print temp_hash
        psm_hash.append(temp_hash)
    return psm_hash


def _get_score_(psm):
    """
    :param psm: peptide to spectrum match dictionairy
    :return: XCorr score (if not available XCorr = 0)
    """
    hit=0
    hit_key=''
    #print psm.keys()
    for key in psm.keys():
        if "xcorr" in key.lower() or 'expectation' in key.lower():
            hit=1
            hit_key=key
    if hit==1:
        return {"XCorr":psm[hit_key]}
    else:
        return {"XCorr":0}

def _get_peptide_sequence_(sequence):
    '''
    :param sequence: peptide sequence with optional tags
    :return: peptide sequence without tags
    '''
    peptide_sequence="".join(re.findall("[a-zA-Z]",sequence))
    return peptide_sequence

def _get_peptide_modifications_(sequence):
    '''
    :param sequence: peptide sequence
    :return: modification parsed out peptide sequence
    '''
    # sequence adjustment for Mascot
    if '|' in sequence:
        split_seq=sequence.split('|')
        split_length=len(split_seq)
        #print split_length
        sequence=''
        i=0
        j=1
        mod_number=0
        mod_value=0
        while i < len(split_seq[0]):
            if mod_number==0:
                inner_split=split_seq[j].split(':')
                mod_number=re.findall(r'\d+', inner_split[0])[0]
                mod_value=inner_split[1]
                #print mod_number
                #print mod_value
            if i < int(mod_number):
                sequence=sequence+split_seq[0][i]
                i+=1
            else:
                sequence=sequence+str(mod_value)
                mod_number=0
                mod_value=0
                if j==split_length-1:
                    mod_number=101
                else:
                    j+=1

    modification=[]
    mod={"position:":"","mass":""}
    hit=0
    pos=0
    for i in range(0,len(sequence)):
        if re.match("[0-9\+\-]",sequence[i]) is not None:
            if hit==0:
                position=pos
                literal_pos=i
                operator=sequence[i]
                value=""
                hit=1
            else:
                value+=sequence[i]
        elif hit==1:
            if sequence[literal_pos-1] in mass.std_aa_mass.keys():
                if operator=="+":
                    value=int(value)+mass.std_aa_mass[sequence[literal_pos-1]]
                else:
                    value=int(value)-mass.std_aa_mass[sequence[literal_pos-1]]
            modification.append({"position":pos,"mass":value})
            hit=0
            pos+=1
        else:
            pos+=1
    return modification

#get_PSM_mzid("/home/vladie/Desktop/Tutorial.test.mzid")
