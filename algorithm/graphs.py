import numpy as np



def Euclidean_distance(x):
    """
    Calculate the distance among each raw of x
    :param x: N X D
                N: the object number
                D: Dimension of the feature
    :return: N X N distance matrix
    """
    x = np.mat(x)
    aa = np.sum(np.multiply(x, x), 1)
    ab = x * x.T
    dist_mat = aa + aa.T - 2 * ab
    dist_mat[dist_mat < 0] = 0
    dist_mat = np.sqrt(dist_mat)
    dist_mat = np.maximum(dist_mat, dist_mat.T)
    dist_mat = np.asarray(dist_mat)
    return dist_mat



def rbf(distance, scale):
    return np.exp(- np.power(distance/scale, 2) )

def computeGraphLaplacian(X, n_neighbors):
    N = X.shape[0]
    distance = Euclidean_distance(X)
    print(distance)
    scale = np.max(distance)
    
    heat_kernel = rbf(distance, scale)
    A = np.zeros([N, N])   #adjacency matrix
    for idx in range(N):
        curr_dist = distance[idx, :]
        curr_sorted = np.argsort(curr_dist)[1:n_neighbors+1]  #dont count itself
        A[idx, curr_sorted] = heat_kernel[idx, curr_sorted]
    A = np.maximum(A, A.transpose())
    D = np.diag(np.sum(A, axis = 0))
    L = D - A
    return L, A, D
