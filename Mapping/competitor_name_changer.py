'''
Created on Mar 2, 2015

@author: bendulken
'''

import os
import sys
import shutil

cell_name=[]
barcode=[]
content=[]
translator=open(sys.argv[1],'rU')
print sys.argv[1]
for line in translator:
    cell_name.append(line.split()[0])
    barcode.append(line.split()[1])
translator.close()

print cell_name
file_name=[]
rootdir=sys.argv[2]
for item in os.listdir(rootdir):
    if os.path.isfile(os.path.join(rootdir, item)):
        file_name.append(os.path.basename(item))

mod_file_name=file_name[1:len(file_name)]
print mod_file_name
new_file_name=[]
for name in mod_file_name:
    file_barcode=name.split("_")[0]
    num=str(name.split("_")[2])
    ind=barcode.index(file_barcode)
    new_name=cell_name[ind]+num
    shutil.copy(rootdir+'/'+name,rootdir+'/new_name/'+new_name)
