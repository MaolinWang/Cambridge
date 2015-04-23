__author__ = 'lin'
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Crippen
from rdkit.Chem.rdMolDescriptors import *
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

file_input = '/home/lin/Downloads/approved.sdf'
file_input_ref = '/home/lin/Downloads/reference.sdf'
file_input_sulfone = '/home/lin/Downloads/sulfone10000.sdf'
file_input_acid = '/home/lin/Downloads/acid10000.sdf'
supplier = Chem.ForwardSDMolSupplier(file_input,sanitize=True)
supplier_ref = Chem.ForwardSDMolSupplier(file_input_ref,sanitize=True)
supplier_sulfone = Chem.ForwardSDMolSupplier(file_input_sulfone,sanitize=True)
supplier_acid = Chem.ForwardSDMolSupplier(file_input_acid,sanitize=True)
remover = SaltRemover()
descriptors=[]
a=0
b=0
c=0
d=0
for x in supplier:
    if x is None:
        continue
    x = remover.StripMol(x,dontRemoveEverything=True)
    descriptors.append([float(Crippen.MolLogP(x)),float(CalcTPSA(x)),float(CalcNumHBA(x)),float(CalcNumHBD(x)),float(CalcExactMolWt(x))])
    a+=1

for x in supplier_ref:
    if x is None:
        continue
    x = remover.StripMol(x,dontRemoveEverything=True)
    descriptors.append([float(Crippen.MolLogP(x)),float(CalcTPSA(x)),float(CalcNumHBA(x)),float(CalcNumHBD(x)),float(CalcExactMolWt(x))])
    b+=1

for x in supplier_sulfone:
    if x is None:
        continue
    x = remover.StripMol(x,dontRemoveEverything=True)
    descriptors.append([float(Crippen.MolLogP(x)),float(CalcTPSA(x)),float(CalcNumHBA(x)),float(CalcNumHBD(x)),float(CalcExactMolWt(x))])
    c+=1

for x in supplier_acid:
    if x is None:
        continue
    x = remover.StripMol(x,dontRemoveEverything=True)
    descriptors.append([float(Crippen.MolLogP(x)),float(CalcTPSA(x)),float(CalcNumHBA(x)),float(CalcNumHBD(x)),float(CalcExactMolWt(x))])
    d+=1

pca = PCA(n_components=2)
X_r = pca.fit(descriptors).transform(descriptors)
print('explained variance ratio (first two components): %s'% str(pca.explained_variance_ratio_))

#plt.figure(figsize=(1024,1024))
for i in range(len(X_r)):
    #if X_r[i][0]<-6000:
        #continue
    if i<a:
        plt.scatter(X_r[i][0], X_r[i][1], label='FDA drugs', color='black')
    if a<=i<(a+b-1):
        plt.scatter(X_r[i][0], X_r[i][1], label='Ref', color='blue')
    if (a+b-1)<=i<(a+b+c-1):
        plt.scatter(X_r[i][0], X_r[i][1], label='sulfone', color='red')
    if (a+b+c-1)<=i<(a+b+c+d-1):
        plt.scatter(X_r[i][0], X_r[i][1], label='acid', color='yellow')
plt.show()
