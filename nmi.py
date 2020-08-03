import numpy as np
from sklearn.metrics.cluster import normalized_mutual_info_score

# merged_ref, merged_own  = np.loadtxt('./test.txt', usecols=range(0,2),unpack=True, skiprows=0) 

merged_ref  = np.loadtxt('./reference/fede_pointMERGED.txt', usecols=range(0,1),unpack=True, skiprows=0) 
merged_own  = np.loadtxt('labels_own.txt', usecols=range(0,1),unpack=True, skiprows=0) 

merged_ref = merged_ref[merged_ref != -9]
merged_own = merged_own[merged_own != -9]

nmi = normalized_mutual_info_score(merged_own, merged_ref)

print(f'NMI: {nmi}')
