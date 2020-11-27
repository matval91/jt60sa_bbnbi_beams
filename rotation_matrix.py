import numpy as np

def rotation_matrix(a,b):
    """
    http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    """
    M = np.zeros((3,3), dtype=float)
    v = np.cross(a,b)
    c = np.dot(a,b)

    vx_mat  = np.array(([0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]))
    vx_mat2 = np.linalg.matrix_power(vx_mat, 2)
    coeff = 1/(1+c)

    M = M+vx_mat+np.multiply(vx_mat2, coeff)+np.identity(3,dtype=float)
    return M
