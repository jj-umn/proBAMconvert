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

from pyteomics import mzid,mass,xml,auxiliary
from fnmatch import fnmatch
import time

#
# Parse mzid file into dictionairy
#

def get_PSM_mzid(psm_file):
    '''
    :param psm_file: mzid file
    :return: dictionairy of parsed mzid file, suitable for proBAMconvert
    '''

    with mzid.read(psm_file) as PSM:
        psm_hash=[]
        accession_hash=_get_accessions_(psm_file)
        for row in PSM:
            temp_hash={"assumed_charge":row['SpectrumIdentificationItem'][0]['chargeState'],"spectrum":row['spectrumID'],"search_hit":[]}
            for psm in row["SpectrumIdentificationItem"]:
                if "Modification" in psm.keys():
                    mod_list=psm['Modification']
                else:
                    mod_list=[]
                proteins=[]
                massdiff=_cal_massdiff_(psm['experimentalMassToCharge'],psm['calculatedMassToCharge'])
                for protein in psm["PeptideEvidenceRef"]:
                        proteins.append({"protein":accession_hash[protein['peptideEvidence_ref']]})
                mod_peptide=_get_mod_peptide_sequence_(psm['peptide_ref'],mod_list)
                modifications=_get_peptide_modifications_(mod_list)
                score=_get_score_(psm)

                temp_hash['search_hit'].append({"hit_rank":psm['rank'],"modifications":modifications,
                                                "modified_peptide":mod_peptide,"peptide":psm['peptide_ref'],
                                                "search_score":{"score":_get_score_(psm),"evalue":_get_evalue_(psm)},
                                                "proteins":proteins,"num_missed_cleavages":"0",
                                                "massdiff":massdiff})
            psm_hash.append(temp_hash)
    return psm_hash

def _get_accessions_(psm_file):
    accession_hash={}
    with open(psm_file,'r') as f:
        count=0
        for line in f:
            count+=1
            if "<PeptideEvidence" in line:
                line=line.split(" ")
                for tag in line:
                    tag=tag.split("=")
                    if tag[0]=="dBSequence_ref":
                        ref=tag[1].replace("\"","")
                    elif tag[0]=="id":
                        id=tag[1].replace("\"","")
                accession_hash[id]=ref
    return accession_hash

def _cal_massdiff_(calc_mass,exp_mass):
    mass_diff=exp_mass-calc_mass
    return mass_diff

def _get_score_(psm):
    """
    :param psm: peptide to spectrum match dictionairy
    :return: XCorr score (if not available XCorr = 0)
    """
    hit=0
    hit_key=''
    score={}
    #print psm.keys()
    for key in psm.keys():
        if "score" in key.lower():
            hit=1
            hit_key=key
    if hit==1:
        return psm[hit_key]
    else:
        return "*"

def _get_evalue_(psm):

    hit=0
    hit_key=''
    score={}
    #print psm.keys()
    for key in psm.keys():
        if "xcorr" in key.lower() or 'expectation' in key.lower() or 'confidence' in key.lower() \
        or "e_value" in key.lower().replace("-","_") or 'evalue' in key.lower():
            hit=1
            hit_key=key
    if hit==1:
        return psm[hit_key]
    else:
        return "*"

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

def get_enzyme_mzid(psm_file):
    f=open(psm_file,'r')
    lines = f.read()
    start = lines.find('<EnzymeName>')
    stop = lines.find('</EnzymeName>')
    line = lines[start:stop]
    line=line.replace('-','_')
    line=line.lower()
    if fnmatch(line,'*trypsin?p*'):
        return 2
    elif fnmatch(line,'*trypsin*'):
        return 1
    elif fnmatch(line, '*lys?c*'):
        return 3
    elif fnmatch(line, '*lys?n*'):
        return 4
    elif fnmatch(line, '*arg?c*'):
        return 5
    elif fnmatch(line, '*asp?n*'):
        return 6
    elif fnmatch(line, '*cnbr*'):
        return 7
    elif fnmatch(line, '*glu?c*'):
        return 8
    elif fnmatch(line, '*pepsina*'):
        return 9
    elif fnmatch(line, '*chymotrypsin*'):
        return 10
    elif fnmatch(line, '*noenzyme*'):
        return 0
    else:
        return "*"

def get_enzyme_specificity_mzid(psm_file):
    f=open(psm_file,'r')
    lines = f.read()
    start = lines.find('<EnzymeName>')
    stop = lines.find('</EnzymeName>')
    line = lines[start:stop]
    line=line.replace('-','_')
    line=line.lower()
    if fnmatch(line,'*semiSpecific=1*'):
        return 1
    elif fnmatch(line,'*semiSpecific=0*'):
        return 2
    else:
        return 6

def extract_comments_from_mzid(psm_file):
    f=open(psm_file,'r')
    comments=[]
    lines = f.read()

    start = lines.find('<AnalysisSoftwareList>')
    stop = lines.find('</AnalysisSoftwareList>')
    line = lines[start:stop]
    for comment in line.split("\n"):
        comments.append(comment)

    start = lines.find('<AnalysisProtocolCollection>')
    stop = lines.find('</AnalysisProtocolCollection>')
    line = lines[start:stop]
    for comment in line.split("\n"):
        comments.append(comment)

    start = lines.find('<AnalysisCollection>')
    stop = lines.find('</AnalysisCollection>')
    line = lines[start:stop]
    for comment in line.split("\n"):
        comments.append(comment)
    return comments

#get_PSM_mzid("/home/vladie/Desktop/proBAMconvert/PeptideShaker.mzid")
# extract_comments_from_mzid("/home/vladie/Desktop/proBAMconvert/Mudpit.mzid")
#get_PSM_mzid("/home/vladie/Desktop/mESC_ignolia/NtermCofr/NtermCofr.mzid")
#get_PSM_mzid("/home/vladie/Desktop/mESC_ignolia/NtermCofr/NtermCofr2.mzid")