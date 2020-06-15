import numpy as np
import matplotlib.pyplot as plt


# load aggregation data set
data_agg = np.loadtxt('./data/aggregation.txt', delimiter='	', usecols=range(0,3), unpack=True)
data_agg_gt = data_agg[-1]
data_agg = data_agg[:-1]
n = data_agg.shape[1]
print(f'data shape: {data_agg.shape}')


clust_id_agg = np.loadtxt('./labels.txt', delimiter='	', usecols=range(0,1), unpack=True)

print(f'labels shape: {clust_id_agg.shape}')

plt.figure()
plt.scatter(data_agg[0], data_agg[1], c = clust_id_agg[0:788])
# plt.scatter(peaks_dp_agg[:,0], peaks_dp_agg[:,1], s = 30, c = 'r')
# plt.scatter(centroids_dp_agg[:,0], centroids_dp_agg[:,1], s = 30, c = 'orange')
plt.tight_layout()
plt.show()