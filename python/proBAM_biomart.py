from __future__ import division

__author__ = 'Volodimir Olexiouk'

from bioservices import BioMart
import sys

#
#Function that links the correct database archive with version number
#
def _get_ensembl_archive_(version):
    '''
    :param version: Ensembl version
    :return: ENSEMBL repository for a specific version
    '''
    version=int(version)
    d={}
    d[83]="www.ensembl.org"
    d[82]="sep2015.archive.ensembl.org"
    d[81]="jul2015.archive.ensembl.org"
    d[80]="may2015.archive.ensembl.org"
    d[79]="mar2015.archive.ensembl.org"
    d[78]="dec2014.archive.ensembl.org"
    d[77]="oct2014.archive.ensembl.org"
    d[76]="aug2014.archive.ensembl.org"
    d[75]="feb2014.archive.ensembl.org"
    d[74]="dec2013.archive.ensembl.org"
    d[73]="sep2013.archive.ensembl.org"
    d[72]="jun2013.archive.ensembl.org"
    d[71]="apr2013.archive.ensembl.org"
    d[70]="jan2013.archive.ensembl.org"
    d[69]="oct2012.archive.ensembl.org"
    d[68]="jul2012.archive.ensembl.org"
    d[67]="may2012.archive.ensembl.org"
    d[66]="feb2012.archive.ensembl.org"
    d[65]="dec2011.archive.ensembl.org"
    d[64]="sep2011.archive.ensembl.org"
    d[63]="jun2011.archive.ensembl.org"
    d[62]="apr2011.archive.ensembl.org"
    d[61]="feb2011.archive.ensembl.org"
    d[60]="nov2010.archive.ensembl.org"
    d[59]="aug2010.archive.ensembl.org"
    d[58]="may2010.archive.ensembl.org"
    d[57]="mar2010.archive.ensembl.org"
    d[56]="sep2009.archive.ensembl.org"
    d[55]="jul2009.archive.ensembl.org"
    d[54]="may2009.archive.ensembl.org"
    if version in d.keys():
        return d[version]
    else:
        raise ValueError('unsupported ensembl version')


#
# Function to map spesies to ENSEMBL dataset names
#
def _get_ensembl_dataset_(species):
    '''
    :param species: full species name
    :return: ensembl species name
    '''
    d={}
    d['homo_sapiens']='hsapiens_gene_ensembl'
    d['mus_musculus']='mmusculus_gene_ensembl'
    d['drosophila_melanogaster']='dmelanogaster_gene_ensembl'
    d['danio_rerio']='drerio_gene_ensembl'

    if species not in d.keys():
        print 'Error: unsupported species'
        print 'Currently supported species:'
        print d.keys()
        raise ValueError('unsupported species')
    return d[species]
#
# Function to create XML readable transcript_id query string
#
def _id_in_xml_query_(transcipt_id):
    '''
    :param transcipt_id: list of transcrip IDs
    :return: XML readable transcript ID query string
    '''
    query=""
    for tr in transcipt_id:
        query+=(str(tr)+",")
    query=query[:-1]
    return query

#
# Function that retrieves cds,strand,chr and ensembl_transcript_id from BioMart
#
def retrieve_data_from_biomart(version,species,transcript_id):
    '''
    :param version: Database version
    :param species: Full species name
    :param transcript_id: list of transcript IDs
    :return: BioMart results
    '''
    tr_query=_id_in_xml_query_(transcript_id)
    version=_get_ensembl_archive_(version)
    dataset=_get_ensembl_dataset_(species)
    biomart = BioMart(host=version)
    biomart.add_dataset_to_xml(dataset)
    biomart.add_filter_to_xml("ensembl_transcript_id",tr_query)
    biomart.add_attribute_to_xml('ensembl_transcript_id')
    biomart.add_attribute_to_xml("chromosome_name")
    biomart.add_attribute_to_xml("strand")
    biomart.add_attribute_to_xml("coding")
    attributes=biomart.attributes(dataset)
    xml_query=biomart.get_xml()
    result=biomart.query(xml_query)
    result=result.split("\n")
    return result

#
# Function that maps Identifiers to ENSEMBl
#
def id_map_ensembl(to_annotation,version,species,psm_protein_id):
    '''
    :param to_annotation: target identifier annotation (i.e. uniprot_swissprot)
    :param version: Database version
    :param species: Full species name
    :param psm_protein_id: list of IDs to be converted
    :return: BioMart results
    '''
    new_psm_protein_id=[]
    query_string=_id_in_xml_query_(psm_protein_id)
    version=_get_ensembl_archive_(version)
    dataset=_get_ensembl_dataset_(species)
    biomart = BioMart(host=version)
    biomart.add_dataset_to_xml(dataset)
    biomart.add_filter_to_xml(to_annotation,query_string)
    biomart.add_attribute_to_xml("ensembl_transcript_id")
    biomart.add_attribute_to_xml("transcript_length")
    biomart.add_attribute_to_xml("uniprot_swissprot")
    attributes=biomart.attributes(dataset)
    xml_query=biomart.get_xml()
    result=biomart.query(xml_query)
    result=result.split("\n")
    return result






