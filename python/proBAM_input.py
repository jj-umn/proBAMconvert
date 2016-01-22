import proBAM_mzid
import proBAM_mzTab
from pyteomics import pepxml
import re


#
# Import PSM file and parse into dictionairy
#
def get_PSM_hash(psm_file,decoy_annotation):
    '''
    :param psm_file: psm file (pepxml,mzid or mztab)
    :raise: IO error: unrecognized file format
    :return: dictionairy of parsed psm file,
    '''
    print 'Importing PSM file'

    # catch mzid file format and parse
    if re.match('^.*\.(mzid)$',psm_file.lower())!=None:
        PSM=proBAM_mzid.get_PSM_mzid(psm_file)

    # catch pepxml file format and parse
    elif re.match('^.*\.(pepxml|pep.xml|xml)$',psm_file.lower())!=None:
        PSM=[]
        PEP=pepxml.read(psm_file,read_schema=False,iterative=True)
        # count = 0
        # parse tags out of protein IDs
        for row in PEP:
            #print row
            #if count==5500:
            #    break
            PSM.append(row)
            #count+=1
        del PEP

    # catch mztab file format and parse
    elif re.match('^.*\.(mztab)$',psm_file.lower())!=None:
        PSM=proBAM_mzTab.get_PSM_mztab(psm_file)

    else:
        raise IOError('Unrecognized file extension, \n ' \
              'Accepted file extensions: .mzid/.pepxml/.pep.xml/.xml')

    # Debug section
    '''
    count=0
    PSM=[]
    for row in PEP:
        print row
        if count==10000:
            break
        count+=1
        PSM.append(row)
    del PEP
    '''
    return PSM
