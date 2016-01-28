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

from pyteomics import mzid,mass,xml


#
# Parse mzid file into dictionairy
#

def get_PSM_mzid(psm_file):
    '''
    :param psm_file: mzid file
    :return: dictionairy of parsed mzid file, suitable for proBAMconvert
    '''
    PSM=mzid.read(psm_file,retrieve_refs=True)
    psm_hash=[]
    count=0

    for row in PSM:
        #print row
        temp_hash={"assumed_charge":row['SpectrumIdentificationItem'][0]['chargeState'],"spectrum":row['spectrumID'],"search_hit":[]}
        for psm in row["SpectrumIdentificationItem"]:
            if "Modification" in psm.keys():
                mod_list=psm['Modification']
            else:
                mod_list=[]
            proteins=[]
            for protein in psm["PeptideEvidenceRef"]:
                    proteins.append({"protein":protein['accession']})
            mod_peptide=_get_mod_peptide_sequence_(psm['PeptideSequence'],mod_list)
            modifications=_get_peptide_modifications_(mod_list)
            score=_get_score_(psm)

            temp_hash['search_hit'].append({"hit_rank":psm['rank'],"modifications":modifications,
                                            "modified_peptide":mod_peptide,"peptide":psm['PeptideSequence'],
                                            "search_score":score,"proteins":proteins,"num_missed_cleavages":"0"})
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

def _get_mod_peptide_sequence_(sequence,modification):
    '''
    :param sequence: peptide sequence with optional tags
    :return: peptide sequence without tags
    '''
    mod_peptide_sequence=""
    for position in range(0,len(sequence)):
        mod_peptide_sequence+=sequence[position]
        for mod in modification:
            if int(mod['location'])==position:
                mod_peptide_sequence+='['+str(mod['monoisotopicMassDelta'])+']'

    return mod_peptide_sequence

def _get_peptide_modifications_(modifications):
    '''
    :param sequence: peptide sequence
    :return: modification parsed out peptide sequence
    '''
    modification=[]
    for mod in modifications:
        pos=mod['location']
        value=mod['monoisotopicMassDelta']
        modification.append({"position":pos,"mass":value})
    return modification

#get_PSM_mzid("/home/vladie/Desktop/mESC_ignolia/NtermCofr/NtermCofr.mzid")
#get_PSM_mzid("/home/vladie/Desktop/mESC_ignolia/NtermCofr/NtermCofr2.mzid")