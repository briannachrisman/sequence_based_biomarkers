import sys
import numpy as np
from scipy.special import comb
from itertools import combinations

order = int(sys.argv[1])
person_taxa_filename = sys.argv[2]
taxa_variant_filename = sys.argv[3]
output_dir = sys.argv[4]

people = []
with open(person_taxa_filename, 'r') as f:
    taxa = next(f).strip().split('\t')
    for line in f:
        people.append(line.strip().split('\t', maxsplit=1)[0])

person_taxa = np.loadtxt(person_taxa_filename, delimiter='\t', skiprows=1, usecols=list(range(1, len(taxa)+1)))

taxa = []
with open(taxa_variant_filename, 'r') as f:
	variants = next(f).strip().split('\t')
	for line in f:
		taxa.append(line.strip().split('\t', maxsplit=1)[0])

taxa_variant = np.loadtxt(taxa_variant_filename, delimiter='\t', skiprows=1, usecols=list(range(1, len(variants)+1)), dtype=bool)
print('Data loaded', 'people', len(people), 'taxa', len(taxa), 'variants', len(variants))


m, n = taxa_variant.shape
cached_combs = np.zeros((n, order+1), dtype=int)
for i in range(n):
	cached_combs[i, :] = [comb(i, x) for x in range(order+1)]

r = int(comb(n, order))
print('num biomarkers', r)

# find sparsity pattern
biomarker_exists = np.zeros((r,), dtype=bool)
for taxa_index in range(len(taxa)):
	var = np.where(taxa_variant[taxa_index, :])[0]
	p = var.shape[0]
	num_indices = cached_combs[p, order]
	if order > 1:
		new_biomarkers = np.fromiter(combinations(var, order), dtype=np.dtype(','.join(['i']*order)), count=int(comb(p, order))).view(np.dtype('i')).reshape(-1, order)
	else:
		new_biomarkers = var[:, np.newaxis]
		
	biomarker_indices = (r-1)*np.ones((num_indices,), dtype=int)
	for j in range(order):
		biomarker_indices -= cached_combs[n-1-new_biomarkers[:, j], order-j]
	biomarker_exists[biomarker_indices] = True

	if taxa_index%100==0:
		print(taxa_index)
print('fraction of biomarkers that exist:', np.sum(biomarker_exists)/r)

biomarker_mapping = np.cumsum(biomarker_exists)-1
np.save('%sbiomarker_exists%d' % (output_dir,order), biomarker_exists)

# now calculate person_biomarker
person_biomarker = np.zeros((len(people), np.sum(biomarker_exists)), dtype=np.float64)
	
for taxa_index in range(len(taxa)):
	#print(taxa_index, end=' ')
	var = np.where(taxa_variant[taxa_index, :])[0]
	p = var.shape[0]
	num_indices = cached_combs[p, order]

	if order > 1:
		new_biomarkers = np.fromiter(combinations(var, order), dtype=np.dtype(','.join(['i']*order)), count=int(comb(p, order))).view(np.dtype('i')).reshape(-1, order)
	else:
		new_biomarkers = var[:, np.newaxis]
		
	biomarker_indices = (r-1)*np.ones((num_indices,), dtype=int)
	for j in range(order):
		biomarker_indices -= cached_combs[n-1-new_biomarkers[:, j], order-j]

	person_indices = np.where(person_taxa[:, taxa_index])[0]
	person_biomarker[np.ix_(person_indices, biomarker_mapping[biomarker_indices])] += np.outer(person_taxa[person_indices, taxa_index], np.ones((num_indices,)))

	if taxa_index%1000==0:
		print(taxa_index)
np.save('%sperson_variant%d_condensed' %  (output_dir,order), person_biomarker)

# write biomarker names
biomarkers = combinations(range(len(variants)), order)
biomarker_exists = np.load(output_dir + 'biomarker_exists%d.npy' % order)
with open('%sbiomarkers%d.txt' % (output_dir,order), 'w+') as f:
    for exists in biomarker_exists:
        biomarker = next(biomarkers)
        if exists:
            f.write('\t'.join([variants[x][1:-1] for x in biomarker]) + '\t' + '\t'.join([str(x) for x in biomarker]) + '\n')

print('%sbiomarkers%d.txt' % (output_dir,order))
print("done! -- success")
