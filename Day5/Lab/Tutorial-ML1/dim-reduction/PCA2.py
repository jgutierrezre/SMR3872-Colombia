import numpy as np
import dask.array as da
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import eigh
from sklearn.datasets import make_moons
import time


def rbf_kernel_pca(X, gamma, n_components):
    """
    RBF kernel PCA implementation.
    Parameters
    ------------
    X: {NumPy ndarray}, shape = [n_examples, n_features]
    gamma: float
    T   uning parameter of the RBF kernel
    n_components: int
        Number of principal components to return
    Returns
    ------------
    X_pc: {NumPy ndarray}, shape = [n_examples, k_features]
        Projected dataset
    """
    # Calculate pairwise squared Euclidean distances
    # in the MxN dimensional dataset.
    sq_dists = pdist(X, "sqeuclidean")
    # Convert pairwise distances into a square matrix.
    mat_sq_dists = squareform(sq_dists)
    # Compute the symmetric kernel matrix.
    K = np.exp(-gamma * mat_sq_dists)
    # Center the kernel matrix.
    N = K.shape[0]
    one_n = np.ones((N, N)) / N
    K = K - one_n.dot(K) - K.dot(one_n) + one_n.dot(K).dot(one_n)
    # Obtaining eigenpairs from the centered kernel matrix
    # scipy.linalg.eigh returns them in ascending order
    eigvals, eigvecs = eigh(K)
    eigvals, eigvecs = eigvals[::-1], eigvecs[:, ::-1]
    # Collect the top k eigenvectors (projected examples)
    X_pc = np.column_stack([eigvecs[:, i] for i in range(n_components)])
    return X_pc


def dask_rbf_kernel_pca(X, gamma, n_components):
    # Convert X to a Dask array if it isn't one already
    if not isinstance(X, da.Array):
        X = da.from_array(X, chunks=(X.shape[0] // 4, X.shape[1]))

    # Calculate square distances
    sq_dists = pdist(X.compute(), "sqeuclidean")
    mat_sq_dists = squareform(sq_dists)

    # Create a Dask array from the square distance matrix
    mat_sq_dists = da.from_array(
        mat_sq_dists, chunks=(mat_sq_dists.shape[0] // 4, mat_sq_dists.shape[1] // 4)
    )

    # Calculate the kernel matrix
    K = da.exp(-gamma * mat_sq_dists)

    # Center the kernel matrix
    N = K.shape[0]
    one_n = da.ones((N, N), chunks=(N // 4, N // 4)) / N
    K = K - one_n.dot(K) - K.dot(one_n) + one_n.dot(K).dot(one_n)

    # Compute eigenvalues and eigenvectors
    eigvals, eigvecs = eigh(K.compute())

    # Sort eigenvalues and eigenvectors in descending order
    eigvals, eigvecs = eigvals[::-1], eigvecs[:, ::-1]

    # Collect the top `n_components` eigenvectors
    X_pc = da.concatenate([eigvecs[:, i : i + 1] for i in range(n_components)], axis=1)

    return X_pc

X, y = make_moons(n_samples=10000, random_state=123)

start_time = time.time()
X_kpca = rbf_kernel_pca(X, gamma=15, n_components=2)
print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()
X_kpca = dask_rbf_kernel_pca(X, gamma=15, n_components=2)
print("--- %s seconds ---" % (time.time() - start_time))
