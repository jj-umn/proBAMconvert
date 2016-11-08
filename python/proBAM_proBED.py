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
from proBAM_coref import *

#######################
### GETTERS/SETTERS ###
#######################
#
# Open file to write to
#
def open_bed_file(directory,name):
    '''
    :param directory: working directory
    :param name: output file name
    :return: return output file IO
    '''
    file=open(directory+name+'.pro.bed',"w")
    return file

#
# Function to write temp_result to file
#
def _write_psm_(temp_result,file):
    '''
    :param temp_result: sam results
    :param file: sam file
    :return: write to sam file IO
    '''
    for tab in temp_result:
        file.write(str(tab)+'\t')
    file.write('\n')

#
# Convert PSM to SAM
#
def PSM2BED(psm_hash,transcript_hash,exon_hash,decoy_annotation,allowed_mismatches,file,map_decoy,rm_duplicates,
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
    print "Commencing generation of BED file"
    # psm_hash.reset()
    if rm_duplicates == "Y":
        dup = {}
    unique_name_generator={}
    for psm in psm_hash:
        # convert unmapped PSMs with their own converter
        if 'search_hit' not in psm.keys():
            continue
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
                                if ("DECOY_"+key) in unique_name_generator:
                                    unique_name_generator["DECOY_"+key] += 1
                                else:
                                    unique_name_generator["DECOY_"+key] = 0
                                temp_result= decoy_PSM_to_BED(psm,row,key,transcript_hash,exon_hash,map_decoy,
                                                              three_frame_translation,
                                                              unique_name_generator["DECOY_"+key],allowed_mismatches)

                                if rm_duplicates=="Y" and temp_result!=False:
                                    print temp_result
                                    dup_key= str(temp_result[0])+"_"+str(temp_result[1])+"_"+str(temp_result[13])
                                    if dup_key not in dup.keys():
                                        dup[dup_key]=1
                                        _write_psm_(temp_result,file)
                                elif temp_result!=False:
                                    _write_psm_(temp_result,file)

                    if decoy==0:
                        key=row['proteins'][p]['protein']
                        if key in unique_name_generator:
                            unique_name_generator[key]+=1
                        else:
                            unique_name_generator[key]=0
                        # Filter out PSM where transcript sequences were not found/ non-existent
                        if key not in transcript_hash.keys():
                            continue
                        # transcript not on an canonical transcript
                        # TODO do this nicer by fetching canonical chr
                        elif len(transcript_hash[key]['chr']) > 4:
                            continue
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
                                continue
                            else:
                                # map peptide on protein and retrieve hit position, iterate over all hits
                                for phit in protein_hit:
                                    temp_result=[None]*25
                                    genome_position=calculate_genome_position(phit[0],
                                                                             transcript_hash[key]['strand'],
                                                                             transcript_hash[key]['5UTR_offset'],
                                                                             transcript_hash[key]['start_exon_rank'],
                                                                             row['peptide'],
                                                                             exon_hash[transcript_hash[key]['transcript_id']],
                                                                             transcript_hash[key]['chr'],
                                                                             three_frame_translation)
                                    CIGAR=compute_cigar(genome_position,
                                                                 exon_hash[transcript_hash[key]['transcript_id']],
                                                                 transcript_hash[key]['strand'],row['peptide'])[0]
                                    #
                                    # Mandatory columns adapted from BED format
                                    #
                                    #chrom
                                    temp_result[0]=str(transcript_hash[key]['chr'])
                                    #chromStart
                                    temp_result[1]=genome_position
                                    #chromStop
                                    temp_result[2] = _calculate_stop_(temp_result[1],CIGAR)
                                    #unique protein accession
                                    temp_result[3]=str(key)+"_"+str(unique_name_generator[key])
                                    #score
                                    temp_result[4]=1000
                                    #strand
                                    if str(transcript_hash[key]['strand'])==str(1):
                                        temp_result[5]="+"
                                    else:
                                        temp_result[5]="-"
                                    #chromstart
                                    temp_result[6]=temp_result[1]
                                    #chromstop
                                    temp_result[7]=temp_result[2]
                                    #reserved
                                    temp_result[8]=0
                                    #blockcount
                                    temp_result[9]=str(len(CIGAR.split('M'))-1)
                                    #list of black sizes
                                    temp_result[10]=_get_block_sizes_(CIGAR)
                                    #list of block starts
                                    temp_result[11]=_get_block_starts_(temp_result[1],CIGAR)
                                    #
                                    #Mandatory proteomics specific columns added to the proBED format
                                    #
                                    #protein accession
                                    temp_result[12]=key
                                    #peptide sequence
                                    temp_result[13]=row['peptide']
                                    #uniqueness
                                    temp_result[14]="."
                                    #genome reference version
                                    temp_result[15]="."
                                    #psm score
                                    temp_result[16]=str(row['search_score']['score'])
                                    #fdr
                                    temp_result[17]=str(row['search_score']['evalue'])
                                    #modifications
                                    temp_result[18]=row['modifications']
                                    #charge
                                    temp_result[19]=psm['assumed_charge']
                                    #exp-mass_to-charge
                                    temp_result[20]="."
                                    #calc-mass-to-charge
                                    temp_result[21]="."
                                    #rank
                                    temp_result[22]=row['hit_rank']
                                    #dataset ID
                                    temp_result[23]=psm['spectrum']
                                    #url
                                    temp_result[24]="."
                                    for i in range(0, len(temp_result)):
                                        if temp_result[i] == '*' or temp_result[i] == []:
                                            temp_result[i] == "."
                                    # remove duplicates if rm_duplicates=Y
                                    if rm_duplicates=="Y":
                                        dup_key= str(temp_result[0])+"_"+\
                                                          str(temp_result[1])+"_"+str(temp_result[13])
                                        if dup_key not in dup.keys():
                                            dup[dup_key]=1
                                            _write_psm_(temp_result,file)
                                    else:
                                        _write_psm_(temp_result,file)
    file.close()

#
# Function to convert decoy PSM to SAM format
#

def decoy_PSM_to_BED(psm, row, key, transcript_hash, exon_hash, map_decoy,three_frame_translation,
                     unique_nr, allowed_mismatches):
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
        return False
    else:
        # map peptide on protein and retrieve hit position, iterate over all hits
        for phit in protein_hit:
            temp_result = [None] * 25
            CIGAR=compute_cigar(temp_result[3],
                          exon_hash[transcript_hash[key]['transcript_id']],
                          transcript_hash[key]['strand'], row['peptide'])[0]
            #
            # Mandatory columns adapted from BED format
            #
            # chrom
            temp_result[0] = str(transcript_hash[key]['chr'])
            # chromStart
            temp_result[1] = calculate_genome_position(phit[0],
                                                     transcript_hash[key]['strand'],
                                                     transcript_hash[key]['5UTR_offset'],
                                                     transcript_hash[key]['start_exon_rank'],
                                                     row['peptide'][::-1],
                                                     exon_hash[transcript_hash[key]['transcript_id']],
                                                     transcript_hash[key]['chr'])
            # chromStop
            temp_result[2] = _calculate_stop_(temp_result[1], CIGAR)
            # unique protein accession
            temp_result[3] = "DECOY_"+str(key) + "_" + str(unique_nr)
            # score
            temp_result[4] = 1000
            # strand
            if str(transcript_hash[key]['strand']) == str(1):
                temp_result[5] = "+"
            else:
                temp_result[5] = "-"
            # chromstart
            temp_result[6] = temp_result[1]
            # chromstop
            temp_result[7] = temp_result[2]
            # reserved
            temp_result[8] = 0
            # blockcount
            temp_result[9] = str(len(CIGAR.split('M'))-1)
            # list of black sizes
            temp_result[10] = _get_block_sizes_(CIGAR)
            # list of block starts
            temp_result[11] = _get_block_starts_(temp_result[1], CIGAR)
            #
            # Mandatory proteomics specific columns added to the proBED format
            #
            # protein accession
            temp_result[12] = key
            # peptide sequence
            temp_result[13] = row['peptide']
            # uniqueness
            temp_result[14] = "."
            # genome reference version
            temp_result[15] = "."
            # psm score
            temp_result[16] = str(row['search_score']['score'])
            # fdr
            temp_result[17] = str(row['search_score']['evalue'])
            # modifications
            temp_result[18] = row['modifications']
            # charge
            temp_result[19] = psm['assumed_charge']
            # exp-mass_to-charge
            temp_result[20] = "."
            # calc-mass-to-charge
            temp_result[21] = "."
            # rank
            temp_result[22] = row['hit_rank']
            # dataset ID
            temp_result[23] = psm['spectrum']
            # url
            temp_result[24] = "."
            for i in range(0,len(temp_result)):
                if temp_result[i]=='*' or temp_result[i]==[]:
                    temp_result[i]=="."
            return temp_result
#
# calculate stop
#
def _calculate_stop_(pos,flag):
    length=0
    number_parsing=''
    for char in flag:
        if char.isdigit()==True:
            number_parsing+=char
        else:
            length+=int(number_parsing)
            number_parsing=''
    return int(pos)+length
#
# get block sized
#
def _get_block_sizes_(flag):
    block_sizes=''
    number_parsing=''
    for char in flag:
        if char.isdigit()==True:
            number_parsing+=char
        elif char=="M":
            if block_sizes=='':
                block_sizes+=number_parsing
            else:
                block_sizes += ','+number_parsing
            number_parsing=''
        else:
            number_parsing=''
    return block_sizes
#
# Get block starts
#
def _get_block_starts_(pos,flag):
    sizes=flag.split('M')
    start_sizes=str(pos)
    length=int(pos)
    number_parsing=''
    if len(sizes)<=2:
        return start_sizes
    else:
        for char in flag:
            if char.isdigit() == True:
                number_parsing += char
            elif char == 'N':
                length += int(number_parsing)
                start_sizes+=','+str(length)
            else:
                length += int(number_parsing)
        return start_sizes


#
# Create SAM header
#
def create_BED_header(file,database,database_v,command_line,psm_file,comments):
    '''
    :param file: output file
    :param version: proBAMconvert version
    :param database: database name
    :param sorting_order: SAM sorting order
    :param database_v: database version
    :return:
    '''
    print "Create BED header"
    header=[]
    header.append('# proBED-version\t1.0')

    header.append('#\t'+str(command_line))
    header.append('#\t database:'+str(database)+'\tdatabaseversion:'+str(database_v))
    # get comments and append comments to file
    if comments!=[]:
        for comment in comments:
            comment=str(comment).rstrip()
            comment=re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', comment, flags=re.M)
            header.append('#'+str(comment))
    comments=extract_comments(psm_file)
    if comments!=[]:
        for comment in comments:
            comment=str(comment).rstrip()
            comment=re.sub(r'(^[ \t]+|[ \t]+(?=:))', '', comment, flags=re.M)
            header.append('#'+str(comment))
    for row in header:
        file.write(row+'\n')

