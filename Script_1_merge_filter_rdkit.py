#/usr/bin/python
__author__ = 'mw648'
from rdkit import Chem
import os

def remover(mol):
    removers=Chem.SaltRemover()
    molnew = removers.StripMol(mol)
    return molnew

def merge_dir(dirname,outname):
    #Notice: outname must not be in the dirname!
    try:
        if os.path.abspath(dirname)==os.path.dirname(outname):
            raise NameError
    except NameError:
        print ('The output file must not be in the directory you gave!\n')
    else:
        supplierout=Chem.SDWriter(outname)
        num=1
        for root, dirs, files in os.walk(dirname): # dirs are ignored
            for name in files:
                filename=os.path.join(root,name)
                supplier=Chem.SDMolSupplier(filename)
                for x in supplier:
                    supplierout.write(x)
                    print ('Complete file %d' % num)
                    num+=1

def filter_SMARTS(inputfile,outputfile,smarts):
    num=1
    try:
        if not os.path.exists(inputfile):
            raise IOError
    except IOError:
        print ('The input file does not exist!\n')
    else:
        num=1
        supplier=Chem.SDMolSupplier(inputfile)
        supplierout=Chem.SDWriter(outputfile)
        patt = Chem.MolFromSmarts(smarts)
        for m in supplier:
            m=remover(m)
            if len(m.GetSubstructMatches(patt))==1:
                supplierout.write(m)
            if num%1000==0:
                print ('Filter %d molecules\n' % num)
            num+=1



if __name__=='__main__':
    filter_SMARTS('/home/mw648/Downloads/export.sdf','/home/mw648/Downloads/acid1.sdf','[S]=O')
#[SX4](=[OX1])(=[OX1])[OX2H,OX1H0-]