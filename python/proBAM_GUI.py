__author__ = 'vladie'

from Tkinter import *
from tkFileDialog import *
import os
import sys
import ScrolledText
import time
import proBAM
import proBAM_input
import proBAM_IDparser
from functools import partial

#
# Update log window periodically
#

class Std_redirector(object):
    def __init__(self,widget):
        self.widget = widget

    def write(self,string):
        self.widget.insert(END,string)
        self.widget.see(END)
        proBam_button.update_idletasks()


def _getProjectName_(tk):
    '''
    :param tk: window
    :return: project name
    '''
    global name
    Label(text='Project name').grid(row=2,column=0)
    v=StringVar()

    e=Entry(tk,textvariable=v)
    e.grid(row=2,column=1)
    name=v

def _openFile_():
    '''
    :return: return location of supplied psm file
    '''
    global psm_file
    psm_file = askopenfilename( filetypes = (("All Files","*.*"),("mzIdentML","*.mzid"),("pepxml","*xml")),title=
                            "Please specify file for proBAM conversion")
    if len(psm_file)<50:
        Label(text=psm_file).grid(row=0,column=1)
    else:
        Label(text=psm_file[:50]+"...").grid(row=0,column=1)
    return psm_file

def _getDirectroy_():
    '''
    :return: returns location of working directory
    '''

    global directory
    directory = askdirectory(title=
                            "Please specify the working directory",initialdir=os.getcwd())
    if len(directory)<50:
        Label(text=directory).grid(row=1,column=1)
    else:
       Label(text=directory[:50]+"...").grid(row=1,column=1)
    directory=str(directory)+'/'
    return str(directory)+'/'

def _getSpecies_(tk):
    '''
    :param tk: window
    :return: selected species
    '''
    global species
    Label(text='Select species').grid(row=3,column=0)
    species= StringVar(tk)
    species.set('homo sapiens')

    OptionMenu(tk,species,'homo sapiens','mus musculus','drosophila melanogaster','danio rerio').grid(row=3,column=1)

def _getDatabase_(tk):
    '''
    :param tk: window
    :return: selected database
    '''
    global database
    Label(text='Select database').grid(row=4,column=0)
    database= StringVar(tk)
    database.set('Ensembl')

    OptionMenu(tk,database,'Ensembl').grid(row=4,column=1)

def _getDatabaseVersion_(tk):
    '''
    :param tk: window
    :return: selected database version
    '''
    global database_v
    Label(text='Select database').grid(row=5,column=0)
    database_v= StringVar(tk)
    database_v.set('83')

    OptionMenu(tk,database_v,'83','82','81','80','79','78','77','76','75','74','73','72','71','70','69','68','67','66','65'
               ,'64','63','62','61','60','59','58','57','56','55','54').grid(row=5,column=1)



def _getAllowedMismatches_(tk):
    '''
    :param tk: window
    :return: selected max mismatches
    '''
    global allowed_mismatches
    Label(text='Allowed mismatches').grid(row=6,column=0)
    allowed_mismatches= StringVar(tk)
    allowed_mismatches.set('0')

    OptionMenu(tk,allowed_mismatches,'0','1','2','3','4','5').grid(row=6,column=1)

#
# Global variable declaration
#

def _get_global_arguments_():
    '''
    :return: global variables
    '''
    global decoy_annotation
    global version
    # can be unknown,unsorted, queryname or coordinate, can be specified by user
    global sorting_order
    global psm_file

    decoy_annotation=['REV_','DECOY_','_REVERSED']
    version='1.0'
    # can be unknown,unsorted, queryname or coordinate, can be specified by user
    sorting_order='unknown'

#
# Print selected variables to log
#

def _print_arguments_():
    print 'directory used:         '+directory
    print 'PSM file:               '+psm_file
    print 'species:                '+species.get().replace(' ','_')
    print 'database:               '+database.get().upper()
    print 'database version:       '+str(int(database_v.get()))
    print 'decoy annotation:       '+str(decoy_annotation)
    print 'allowed mismatches:     '+str(int(allowed_mismatches.get()))
    print 'proBAM convert version: '+str(version)
    print 'sorting order:          '+sorting_order
    print 'project name:           '+str(name.get())

#
# Execute proBAMconvert
#
def execute_proBAM(root):
    '''
    :param root: window root
    :return:
    '''
    root.config(cursor="watch")
    root.update()
    start_time = time.time()                                 # start timing function

    # get and print arguments
    _get_global_arguments_()
    _print_arguments_()
    print '***************************************************************************'

    # create connection to SAM file
    file=proBAM.open_sam_file(directory,name.get())

    # hash PSM_DATA and define variables
    psm_hash=proBAM_input.get_PSM_hash(psm_file,decoy_annotation)

    parse_results=proBAM_IDparser.parseID(psm_hash,species.get().replace(' ','_'),
                                       database.get().upper(),decoy_annotation,int(database_v.get()))
    annotation=parse_results[1]
    psm_hash=parse_results[0]
    transcript_hash=annotation[0]
    exon_hash=annotation[1]

    # convert to SAM
    proBAM.create_SAM_header(file,version,database.get().upper(),sorting_order,
                             int(database_v.get()),species.get().replace(' ','_'))
    proBAM.PSM2SAM(psm_hash,transcript_hash,exon_hash,decoy_annotation,int(allowed_mismatches.get()),file)
    proBAM.sam_2_bam(directory,name.get())

    print("proBAM conversion succesful")
    print("%f seconds" % (time.time() - start_time))         # output script run time
    root.config(cursor="")
#
# GUI creation
#

def GUI():
    '''
    :return:
    '''
    #
    # Window initiation
    #
    #os.chdir('/home/vladie/Desktop/proBAM_mzTab')
    root = Tk()
    root.title("proBAM converter: convert pepxml/mzid to proBAM")
    root.geometry("600x700")

    #
    # Rederict standart output to console
    #
    Label(text='Console:').grid(row=10,columnspan=2)
    global std_label
    std_frame=Frame().grid(row=11,columnspan=2)

    std_text = ScrolledText.ScrolledText(std_frame)
    std_text.grid(row=12,columnspan=2)

    sys.stdout = Std_redirector(std_text)
    sys.stderr = Std_redirector(std_text)
    sys.stdin  = Std_redirector(std_text)

    #
    # Create grid for placement
    #
    root.grid()

    #
    # Button and utility assignment
    #

    Button(text='Choose file',command=_openFile_,fg='blue').grid(row=0,column=0)
    Button(text='working directory',command=_getDirectroy_,fg='blue').grid(row=1,column=0)
    _getProjectName_(root)
    species=_getSpecies_(root)
    _getDatabase_(root)
    _getDatabaseVersion_(root)
    _getAllowedMismatches_(root)
    execute_proBAM_argumented=partial(execute_proBAM,root)
    global proBam_button
    proBam_button=Button(text='Convert',fg="blue",command=execute_proBAM_argumented)
    proBam_button.grid(row=7,columnspan=2)
    root.update_idletasks()
    root.mainloop()




 ####################
### MAIN PROGRAM ###
####################

if __name__=='__main__':
    #start GUI
    os.chdir("/home/vladie/Desktop/mESC_ignolia")
    GUI()



