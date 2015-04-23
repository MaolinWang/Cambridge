from rdkit import Chem
from rdkit.Chem import Crippen
import statsmodels.api as sm
import matplotlib.pyplot as plt

file_input = '/home/lin/Downloads/reference.sdf'
file_input_sulfone = '/home/lin/Downloads/sulfone10000.sdf'
file_input_acid = '/home/lin/Downloads/acid10000.sdf'
supplier = Chem.SDMolSupplier(file_input)
supplier_sulfone = Chem.SDMolSupplier(file_input_sulfone)
supplier_acid = Chem.SDMolSupplier(file_input_acid)
mollogp = []
molaf = []
molpre = []
for x in supplier:
    mollogp.append(float(Crippen.MolLogP(x)))
    molpre.append(float(x.GetProp('pre')))
    molaf.append(float(x.GetProp('af')))

X = []
X_sulfone = []
X_acid = []
for x in supplier_sulfone:
    X_acid.append([float(x.GetProp('pre')), float(Crippen.MolLogP(x))])
for x in supplier_acid:
    X_acid.append([float(x.GetProp('pre')), float(Crippen.MolLogP(x))])
Y = []
for i in range(len(mollogp)):
    X.append([molpre[i], mollogp[i]])
    Y.append(molaf[i])
X = sm.add_constant(X)
X_acid = sm.add_constant(X_acid)
model = sm.OLS(Y, X)
results = model.fit()
print(results.summary())
plt.scatter(results.predict(X), Y, color='black',label='Ref',s=80)
plt.scatter(results.predict(X_acid[0]), results.predict(X_acid[0]), color='red',label='sulfone',s=80)
plt.scatter(results.predict(X_acid[1:]), results.predict(X_acid[1:]), color='blue',label='acid',s=80)
plt.legend()
#plt.plot(clf.predict(X), clf.predict(X), color='blue',linewidth=3)
plt.show()