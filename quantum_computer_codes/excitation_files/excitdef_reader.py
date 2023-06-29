import sys,os

working_directory = os.getcwd()
sys.path.insert(0,working_directory)

from parameters import excit_document

"""This function reads a excitation.def file and converts it into a easily manipulated object.
In this case, said object is a list. Each items of this list are the individual lines. 
These items are also themselves lists, each containing in this same order: t ri ra rb."""
def excitdef_reader(document_name,file_location=''):
    file = open(os.path.join(file_location,document_name)).read()
    lines_doc = file.split('\n')[5:-1]
    for line_number in range(len(lines_doc)):
        lines_doc[line_number] = lines_doc[line_number].split() 
        lines_doc[line_number] = [eval(number) for number in lines_doc[line_number]]
    
    return lines_doc
