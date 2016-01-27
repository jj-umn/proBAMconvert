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

import csv
from pyteomics import mass,xml

#
# Parse mztab file into dictionairy
#

def get_PSM_mztab(psm_file):
    '''
    :param psm_file: mztab file
    :return: dictionairy of converted mztab file, suitable for proBAM conversion
    '''
    # connect to modification datbases and parse in desired format (UNIMOD & PSIMOD)
    unimod=_unimod_parser_()
    psimod=_psimod_xml_parser_()
    with open(psm_file,'rU') as csvfile:
        # skip comments
        for row in csvfile:
            tag= row[0:3]
            if tag=='COM' or tag=='MTD' or tag=="" or tag=="\n":
                next(csvfile)
            else:
                break

        #skip header
        next(csvfile)

        #read mztab PSM rows
        mztab=csv.reader(csvfile,delimiter='\t')

        #hash to store variables in processable format
        psm_hash=[]

        #iterate over all iterations and bundle psm from same spectrum
        spectrum={}
        for row in mztab:
            if row!="" and row!="\n" and len(row)!=0:
                if row[0]=="PSM":
                    print row[0]
                    if row[10] not in spectrum.keys():
                        spectrum[row[10]]=[]
                        spectrum[row[10]].append(row)
                    else:
                        spectrum[row[10]].append(row)


        #iterate over all spectrum to store in processable dictionairy
        for key in spectrum.keys():
            if spectrum[key][0][12]!="":
                temp_hash={"assumed_charge":spectrum[key][0][12],"spectrum":key,"search_hit":[]}
            else:
                temp_hash={"assumed_charge":0,"spectrum":key,"search_hit":[]}
            for psm in spectrum[key]:
                proteins=[]
                proteins.append({"protein":psm[3],'peptide_prev_aa':psm[15],"peptide_next_aa":psm[16]})
                modifications=_get_modifications_(psm[1],psm[9],unimod,psimod)
                temp_hash['search_hit'].append({"hit_rank":1,"modifications":modifications,
                                                "modified_peptide":psm[1],"peptide":psm[1],
                                                "search_score":{"XCorr":psm[8]},
                                                "proteins":proteins,"num_missed_cleavages":"0"})
            psm_hash.append(temp_hash)
    return psm_hash
#
# connect to unimod DB, parse and store in dict
#
def _unimod_parser_():
    '''
    :return: dictionairy with unimod IDs as keys and avg mass as value
    '''
    # create UNIMOD dictionairy (record_id->avg mass)
    unimod_dict={}
    # connect UNIMOD database
    unimod_db=mass.Unimod(source='http://www.unimod.org/xml/unimod.xml')
    for mod in  unimod_db.mods:
        unimod_dict[(mod['record_id'])]=mod['avge_mass']
    return unimod_dict

#
# Parse modifications
#
def _get_modifications_(peptide,mods,unimod,psimod):
    '''
    :param peptide: peptide sequence
    :param mods: peptide modification in unimod or psimod format
    :param unimod: unimod dictionairy
    :param psimod: psimod dictionairy
    :return: list of modification dictionairies for this peptide
    '''
    # seperate modification string into list of modifications
    mod_list=mods.split(',')

    # create variable to store results
    modification=[]

    #iterate over all modifications
    for mod in mod_list:
        if mod=="0" or mod.upper()=='NULL':
            break

        #partition unimod variable, separating position from unimod ID
        mod_partitions=mod.split("-")

        #calculate mass
        db_type=mod_partitions[1].split(':')[0]
        if db_type=="UNIMOD":
            #look up avg mass in unimod dict and store modification array
            mass=unimod[int(mod_partitions[1].split(':')[1])]
            modification.append({"position":mod_partitions[0],"mass":mass})
        elif db_type=="MOD":
            mass=psimod[mod_partitions[1]]
            modification.append({"position":mod_partitions[0],"mass":mass})

    return modification

#
# Connect to psimod, parse IDs and store in dict
#
def _psimod_xml_parser_():
    '''
    :return: dictionairy with psimod IDs as key and avg mass as value
    '''
    # create dictionairy to store id/avg difference mass
    psimod={}
    # parse and store psimod data
    import urllib2
    from lxml import etree
    data=urllib2.urlopen("http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/mod/data/PSI-MOD.obo.xml")
    root = etree.parse(data)
    root2 = etree.tostring(root)
    for term in root.xpath('//term'):
        for id in term.xpath('./id'):
            id = id.text
        for xref in term.xpath('./xref_analog'):
            for dbname in xref.xpath('./dbname'):
                if dbname.text=="DiffAvg":
                    avg_diff_mass=xref.xpath('./name/text()')
                    psimod[id]=avg_diff_mass[0]
    return psimod

