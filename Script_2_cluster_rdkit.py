#/usr/bin/python
__author__ = 'Maolin Wang'
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Cluster import ClusterUtils
from rdkit.ML.Cluster import Murtagh

#File to be specified
file_input='/home/lin/work/test/acid_1000.sdf'
file_output='/home/lin/work/test/acid_1000.sdf'

supplier=Chem.SDMolSupplier(file_input)
fps = [FingerprintMols.FingerprintMol(x) for x in supplier]
#Euclidean distance
dists=[]
for i in range(len(fps)):
    for j in range(i):
        dists.append(1.-DataStructs.TanimotoSimilarity(fps[i],fps[j]))
#cluster
clusts = Murtagh.ClusterData(dists,len(fps),Murtagh.WARDS,isDistData=1)
splitClusts=ClusterUtils.SplitIntoNClusters(clusts[0],2000)
centroids =[ClusterUtils.FindClusterCentroidFromDists(x,dists) for x in splitClusts]
#generate new sdf
supplier=Chem.SDMolSupplier(file_input)
num=0
w = Chem.SDWriter(file_output)
for x in supplier:
    if num in centroids:
        w.write(x)
    num+=1


