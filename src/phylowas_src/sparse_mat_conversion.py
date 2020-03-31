from scipy import io
from scipy import sparse
import numpy as np
import pandas as pd
mat = np.load("/oak/stanford/groups/dpwall/users/kpaskov/Biomarkers/person_variant3_condensed.npy")
print("1")
ids = np.random.choice(list(range(np.shape(mat)[1])), 100000, replace=False)
print("2")
mat = mat[:,ids]
names = pd.read_csv("/oak/stanford/groups/dpwall/users/kpaskov/Biomarkers/biomarkers3.txt", header=None, sep='\t')
names = names.iloc[ids,]
print("4")
_, idx = np.unique(mat, axis=1, return_index=True)
print("4")
mat = mat[:,idx]
names = names.iloc[idx,:]
print("HERE")
np.save('/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/yogurt/person_biomarker/DECIPHER/person_variant3_condensed.npy', np.unique(mat, axis=1).astype(np.float64))
names.to_csv('/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/yogurt/person_biomarker/biomarkers3.txt', header=None, index=None)