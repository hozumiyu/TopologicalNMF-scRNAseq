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


def heatKernel(distance, n_neighbors, scale):
    N = distance.shape[0]
    heat_kernel = np.zeros([N, N])
    for idx in range(N):
        curr_dist = distance[idx, :]
        curr_sorted = np.argsort(curr_dist)[1:n_neighbors+1]  #dont count itself
        heat_kernel[idx, curr_sorted] = np.exp(-curr_dist[curr_sorted]**2/ scale)
    return heat_kernel
    
    

def computePersistentLaplacian(X, n_neighbors, weights):
    N = X.shape[0]
    T = len(weights)
    weights = weights / np.sum(weights)
    distance = Euclidean_distance(X)
    scale = np.max(distance)**1.5
    
    HEAT = heatKernel(distance, n_neighbors, scale)
    l_min = np.min(HEAT[HEAT > 0])
    print(l_min)
    l_max = np.max(HEAT)
    d = l_max - l_min
    A = np.zeros([N, N])   #adjacency matrix
    for idx in range(1, T+1):
        A[HEAT >= (idx/T)*d + l_min] += weights[T - idx]
    A = A + A.T - A*A.T
    #A = A + A.T - A*A.T   #symetrize
    D = np.diag(np.sum(A, axis = 0))
    L = D - A
    return L, A, D
