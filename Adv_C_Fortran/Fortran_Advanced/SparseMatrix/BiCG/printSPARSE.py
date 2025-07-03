import numpy as np
from scipy.sparse import csr_matrix

# Original dense matrix in row-major format
dense_data = [
    3., 0., 1., 0., 0.,
    0., 4., 0., 0., 0.,
    0., 7., 5., 9., 0.,
    0., 0., 0., 0., 2.,
    0., 0., 0., 6., 5.
]

n = 5
a = np.array(dense_data).reshape((n, n))

# Threshold for including values (like in sprsin)
thresh = 1.0

# Extract non-zero entries above the threshold
row_idx, col_idx = np.nonzero(np.abs(a) >= thresh)
data = a[row_idx, col_idx]

# Create CSR sparse matrix
sparse_matrix = csr_matrix((data, (row_idx, col_idx)), shape=(n, n))

# Print COO-style (row, col) value
coo = sparse_matrix.tocoo()
for i, j, v in zip(coo.row, coo.col, coo.data):
    print(f"({i+1},{j+1}) {v:.2f}")

print(sparse_matrix)