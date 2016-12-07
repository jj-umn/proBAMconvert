##############################################################################
# Copyright 2016  Volodimir Olexiouk, Gerben Menschaert                      #
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

from Tkinter import *
import ttk
from tkFileDialog import *
import os
import sys
import ScrolledText
import time
import proBAM
import proBAM_input
import proBAM_IDparser
import proBAM_proBED
from functools import partial
import webbrowser
import thread

#import needed for stand-alone executable
import proBAM_mzTab
import proBAM_mzid
import proBAM_ENSEMBL
import proBAM_biomart
import six
import packaging
import packaging.version
import packaging.specifiers
import packaging.markers
import csv
from pyteomics import mass,xml
from pyteomics import mzid,mass,xml
import re
import time
from itertools import imap
import operator
import re
import pysam
import pysam.ctabixproxies
from bioservices import BioMart
import sys

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
    Label(text='Project name',background="#f2f2f2",width=30,anchor=W).grid(row=2,column=0)
    v=StringVar()

    e=Entry(tk,textvariable=v)
    e.grid(row=2,column=1)
    name=v

def _openFile_():
    '''
    :return: return location of supplied psm file
    '''
    global psm_file
    psm_file = askopenfilename(title="Please specify file for proBAM conversion",)
    if len(psm_file)<20:
        Label(text=psm_file,background="#f2f2f2",width=25).grid(row=0,column=1)
    else:
        Label(text='...'+psm_file[-20:],background="#f2f2f2",width=25).grid(row=0,column=1)
    return psm_file

def _getDirectroy_():
    '''
    :return: returns location of working directory
    '''

    global directory
    directory = askdirectory(title=
                            "Please specify the working directory",initialdir=os.getcwd())
    if len(directory)<20:
        Label(text=directory,background="#f2f2f2",width=25).grid(row=1,column=1)
    else:
       Label(text="..."+directory[-20:],background="#f2f2f2",width=25).grid(row=1,column=1)
    directory=str(directory)+'/'
    return str(directory)+'/'

def _getSpecies_(tk):
    '''
    :param tk: window
    :return: selected species
    '''
    global species
    Label(text='Select species',background="#f2f2f2",width=30,anchor=W).grid(row=3,column=0)
    species= StringVar(tk)
    species.set('homo sapiens')

    menu=OptionMenu(tk,species,'homo sapiens','mus musculus','drosophila melanogaster','danio rerio','arabidopsis thaliana')
    menu.config(width=15)
    menu.grid(row=3,column=1)

def _getDatabase_(tk):
    '''
    :param tk: window
    :return: selected database
    '''
    global database
    Label(text='Select database',background="#f2f2f2",width=30,anchor=W).grid(row=4,column=0)
    database= StringVar(tk)
    database.set('Ensembl')

    menu=OptionMenu(tk,database,'Ensembl')
    menu.config(width=15)
    menu.grid(row=4,column=1)

def _getDatabaseVersion_(tk):
    '''
    :param tk: window
    :return: selected database version
    '''
    global database_v
    Label(text='Select database version',background="#f2f2f2",width=30,anchor=W).grid(row=5,column=0)
    database_v= StringVar(tk)
    database_v.set('85')

    menu=OptionMenu(tk,database_v,'85','84','83','82','81','80','79','78','77','76',
                                    '75','74','73','72','71','70','69','68','67','66','65',
                                    '64','63','62','61','60','59','58','57','56','55','54')
    menu.config(width=15,background="#f2f2f2")
    menu.grid(row=5,column=1)

def _getMapDecoy_(tk):
    '''
    :param tk: window
    :return: selected map decoy option
    '''
    global allow_decoys
    Label(text="Allow Decoys",background="#f2f2f2",width=30,anchor=W).grid(row=6,column=0)
    allow_decoys=StringVar(tk)
    allow_decoys.set('N')
    menu=OptionMenu(tk,allow_decoys,'N','Y')
    menu.config(width=15)
    menu.grid(row=6,column=1)

def _getRMDuplicates_(tk):
    '''
    :param tk: window
    :return: selected map duplicates option
    '''
    global rm_duplicates
    Label(text="remove duplicate psm mappings",background="#f2f2f2",width=30,anchor=W).grid(row=7,column=0)
    rm_duplicates=StringVar(tk)
    rm_duplicates.set("N")
    menu=OptionMenu(tk,rm_duplicates,"N","Y")
    menu.config(width=15)
    menu.grid(row=7,column=1)

def _getAllowedMismatches_(tk):
    '''
    :param tk: window
    :return: selected max mismatches
    '''
    global allowed_mismatches
    Label(text='Allowed mismatches',background="#f2f2f2",width=30,anchor=W).grid(row=9,column=0)
    allowed_mismatches= StringVar(tk)
    allowed_mismatches.set('0')

    menu=OptionMenu(tk,allowed_mismatches,'0','1','2','3','4','5')
    menu.config(width=15)
    menu.grid(row=9,column=1)
#
# Center the toplevel widget
#
def center(toplevel):
    toplevel.update_idletasks()
    w = toplevel.winfo_screenwidth()
    h = toplevel.winfo_screenheight()
    size = tuple(int(_) for _ in toplevel.geometry().split('+')[0].split('x'))
    x = w/2 - size[0]/2
    y = h/2 - size[1]/2
    toplevel.geometry("%dx%d+%d+%d" % (size + (x, y)))

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
    global comments
    global pre_picked_annotation
    global three_frame_translation
    global include_unmapped

    comments=[]
    decoy_annotation=['REV_','DECOY_','_REVERSED','REVERSED_','_DECOY']
    version='1.0'
    include_unmapped='Y'
    three_frame_translation="N"
    # can be unknown,unsorted, queryname or coordinate, can be specified by user
    sorting_order='unknown'
    probed='N'

    pre_picked_annotation="First"

def _advanced_settings_(tk):
    global advanced_settings
    execute_open_advanced_settings = partial(execute_proBAM, tk)
    advanced_settings=Button(text='Advanced Settings',fg="black",command=_open_advanced_settings_,width=30,height=1,
                             compound=LEFT,padx=0)
    advanced_settings.grid(row=10,column=0)

def _manual_(tk):

    advanced_settings=Button(text='Manual',fg="black",command=_open_manual_,width=30,height=1,
                             compound=LEFT,padx=0)
    advanced_settings.grid(row=10,column=1)
def _open_manual_():
    import webbrowser
    url="http://probam.biobix.be/manual"
    webbrowser.open(url, new=0, autoraise=True)

def _sortingOrder_(tk):
    '''
    :param tk: window
    :return: selected max mismatches
    '''
    global new_sorting_order
    Label(tk,text='sorting order',background="#f2f2f2",width=30,anchor=W).grid(row=0,column=0)
    new_sorting_order= StringVar(tk)
    new_sorting_order.set('unknown')

    menu=OptionMenu(tk,new_sorting_order,'unknown','unsorted', 'queryname','coordinate')
    menu.config(width=15)
    menu.grid(row=0,column=1)

def _convert_to_probed_(tk):
    global probed
    Label(tk, text='convert to proBED', background="#f2f2f2", width=30, anchor=W).grid(row=8, column=0)
    probed = StringVar(tk)
    probed.set('N')

    menu = OptionMenu(tk, probed, 'Y','N')
    menu.config(width=15)
    menu.grid(row=8, column=1)

def _get3frame_(tk):
    '''
    :param tk: window
    :return: selected map duplicates option
    '''
    global new_three_frame_translation
    Label(tk,text="3-frame translation",background="#f2f2f2",width=30,anchor=W).grid(row=2,column=0)
    new_three_frame_translation=StringVar(tk)
    new_three_frame_translation.set("N")
    menu=OptionMenu(tk,new_three_frame_translation,"N","Y")
    menu.config(width=15)
    menu.grid(row=2,column=1)

def _decoyAnnotation_(tk):
    global new_decoy_annotation
    Label(tk, text='decoy annotation(s)', background="#f2f2f2", width=30, anchor=W).grid(row=1, column=0)
    new_decoy_annotation= StringVar(tk)
    entry=Entry(tk, textvariable=new_decoy_annotation)
    new_decoy_annotation.set('REV_,DECOY_,_REVERSED,REVERSED_,_DECOY')
    entry.grid(row=1,column=1)

def _comments_(tk):
    global new_comments
    Label(tk, text='add comment(s):', pady=5, background="#f2f2f2", width=30, anchor=W).grid(row=5,column=0)
    new_comments = StringVar(tk)
    text=Text(tk)
    text.config(background="white",height=5,width=60)
    text.grid(row=6,columnspan=2)

def _include_unmapped_(tk):
    global new_include_unmapped
    Label(tk, text="include_unmapped", background="#f2f2f2", width=30, anchor=W).grid(row=4, column=0)
    new_include_unmapped = StringVar(tk)
    new_include_unmapped.set("Y")
    menu = OptionMenu(tk, new_include_unmapped, "Y", "N")
    menu.config(width=15)
    menu.grid(row=4, column=1)

def _pre_picked_annotation(tk):
    global new_pre_picked_annotation
    Label(tk,text='annotation identifiers',background="#f2f2f2",width=30,anchor=W).grid(row=3,column=0)
    new_pre_picked_annotation= StringVar(tk)
    new_pre_picked_annotation.set('First')

    menu=OptionMenu(tk,new_pre_picked_annotation,'First','Ensembl_tr', 'Ensembl_pr','UniProt_ACC','UniProt_Entry')
    menu.config(width=15)
    menu.grid(row=3,column=1)

def _save_and_exit_(top):
    global sorting_order
    global decoy_annotation
    global comments
    global pre_picked_annotation
    global include_unmapped
    global three_frame_translation
    if new_three_frame_translation.get()=='Y':
        three_frame_translation='Y'
    if new_sorting_order.get()!='':
        sorting_order=new_sorting_order.get()
    if new_decoy_annotation.get() != '':
        decoy_annotation=new_decoy_annotation.get().split(',')
    if new_comments.get()!='':
        comments=new_comments.get().split("\n")
    if new_pre_picked_annotation.get()!='':
        pre_picked_annotation=new_pre_picked_annotation.get()
    if new_include_unmapped.get()=='N':
        include_unmapped='N'
    top.destroy()

def _open_advanced_settings_():

    #create Toplevel window
    top = Toplevel()
    center(top)
    top.title("Advanced Settings")
    top.configure(background="#f2f2f2", borderwidth=20)
    top.grid()

    # create widgets
    _sortingOrder_(top)
    _decoyAnnotation_(top)
    _include_unmapped_(top)
    _comments_(top)
    _get3frame_(top)
    _pre_picked_annotation(top)

    # create partial save and exit for tk
    save_and_exit_argumented = partial(_save_and_exit_, top)

    save_and_exit_button=Button(top,text='Save & Exit',fg="black",command=save_and_exit_argumented,width=20,height=2,borderwidth=3)
    save_and_exit_button.grid(row=8,columnspan=2,pady=10)

    top.geometry("500x500")


#
# Print selected variables to log
#

def _print_arguments_():
    print 'directory used:          '+directory
    print 'PSM file:                '+psm_file
    print 'species:                 '+species.get().replace(' ','_')
    print 'database:                '+database.get().upper()
    print 'database version:        '+str(int(database_v.get()))
    print 'decoy annotation:        '+str(decoy_annotation)
    print 'allowed mismatches:      '+str(int(allowed_mismatches.get()))
    print 'proBAMconvert version:   '+str(version)
    print 'sorting order:           '+ sorting_order
    print 'project name:            '+str(name.get())
    print 'remove duplicate PSMs:   '+rm_duplicates.get()
    print '3-frame translation:     '+three_frame_translation
    print 'convert to proBED        '+probed.get()
    print 'pre picked annotation    '+pre_picked_annotation
    print 'include unmapped PSMs    '+include_unmapped

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

    start_time = time.time()                                # start timing function

    Label(text='console:', pady=5, width=70, background="#f2f2f2").grid(row=12, columnspan=2)
    std_text = ScrolledText.ScrolledText(root,height=20)
    std_text.grid(row=15,columnspan=2)

    sys.stdout = Std_redirector(std_text)
    sys.stderr = Std_redirector(std_text)
    sys.stdin  = Std_redirector(std_text)

    root.update()
    # get and print arguments
    try:
        _print_arguments_()
        command_line = "python proBAM.py --name " + str(name.get()) + " --mismatches " + str(
            allowed_mismatches.get()) + " --version " + str(database_v.get()) \
                       + " --database " + str(database.get().upper()) + " --species " + str(species.get()) + " --file " + str(psm_file) + \
                       " --directory " + str(directory) + " --rm_duplicates " + str(rm_duplicates.get()) + \
                       " --tri_frame_translation " + \
                       str(three_frame_translation+" --pre_picked_annotation "+pre_picked_annotation) +\
                       " --include_unmapped "+str(include_unmapped)
        print '\n'


        # hash PSM_DATA and define variables
        psm_hash=proBAM_input.get_PSM_hash(psm_file,decoy_annotation)
        parse_results=proBAM_IDparser.parseID(psm_hash,species.get().replace(' ','_'),
                                           database.get().upper(),decoy_annotation,int(database_v.get()),
                                            three_frame_translation,pre_picked_annotation)

        annotation = parse_results[1]
        psm_hash = parse_results[0]
        transcript_hash = annotation[0]
        exon_hash = annotation[1]
        id_map = parse_results[2]

        # convert to SAM
        if probed.get()!='Y':
            file = proBAM.open_sam_file(directory, name.get())
            proBAM.create_SAM_header(file, version, database.get().upper(), sorting_order, database_v.get(),
                                     species.get(), command_line, psm_file,
                              comments)
            proBAM.PSM2SAM(psm_hash, transcript_hash, exon_hash, decoy_annotation, int(allowed_mismatches.get()),
                           file, rm_duplicates.get(),three_frame_translation,psm_file,id_map,root)
            proBAM.compute_NH_XL(directory, name.get(),include_unmapped)
            proBAM.sam_2_bam(directory, name.get())
        else:
            file = proBAM_proBED.open_bed_file(directory, name.get())
            proBAM_proBED.create_BED_header(file, database.get().upper(), database_v.get(), command_line,
                                            psm_file, comments)
            proBAM_proBED.PSM2BED(psm_hash, transcript_hash, exon_hash, decoy_annotation,
                                  int(allowed_mismatches.get()), file, rm_duplicates.get(),
                                  three_frame_translation, id_map, root, database_v.get(), species.get())

        root.config(cursor="")
        print("proBAM conversion succesful")
        print("%f seconds" % (time.time() - start_time))         # output script run time
    except Exception,e:
        print "ERROR:"+str(e)+"\n" \
              "Please check if all parameters were supplied correctly, if the error keeps occuring contact " \
              "the developers at \"https://github.com/Biobix/proBAMconvert/issues\" " \
              "and provide the file to be processed along with the following error message:"
        print e
        root.config(cursor="")
#
# get path of script
#
def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

#
# GUI creation
#

def GUI():
    #
    # Window initiation
    #
    #todo refresh window periodically (now its stuck between processes)
    root = Tk()
    center(root)
    root.title("proBAMconvert")
    root.geometry("620x725")
    root.configure(background="#f2f2f2",borderwidth=20)
    logo=PhotoImage(file=getScriptPath()+'/'+"proBAMconvert_logo.gif")
    Label(root, image=logo, background="#f2f2f2").grid(row=15,columnspan=2)

    #style
    import tkFont
    default_font = tkFont.nametofont("TkDefaultFont")
    default_font.configure(size=9,weight=tkFont.BOLD,family="MS Free Sans")
    std_frame=Frame().grid(row=13,columnspan=2)

    # get scrolled text
    Label(text='', pady=5, width=70, background="#f2f2f2").grid(row=12, columnspan=2)
    global std_label

    #
    # Create grid for placement
    #
    root.grid()

    #
    # Button and utility assignment
    #

    Button(text='Choose file',command=_openFile_,fg='#0099cc',width=30,anchor=W,justify=CENTER).grid(row=0,column=0)
    Button(text='working directory',command=_getDirectroy_,fg='#0099cc',width=30,anchor=W,
           justify=CENTER).grid(row=1,column=0)
    _getProjectName_(root)
    species=_getSpecies_(root)
    _getDatabase_(root)
    _getDatabaseVersion_(root)
    _getRMDuplicates_(root)
    _convert_to_probed_(root)
    _getAllowedMismatches_(root)
    _advanced_settings_(root)
    _manual_(root)
    execute_proBAM_argumented=partial(execute_proBAM,root)
    global proBam_button
    proBam_button=Button(text='Convert',fg="#0099cc",command=lambda:
        thread.start_new_thread(execute_proBAM_argumented, ()),width=20,height=2,borderwidth=3)
    proBam_button.grid(row=11,columnspan=2,pady=10)
    root.update_idletasks()
    root.mainloop()


####################
### MAIN PROGRAM ###
####################

if __name__=='__main__':
    #start GUI
    os.chdir("/home/vladie/Desktop/proBAMconvert")
    _get_global_arguments_()
    GUI()



