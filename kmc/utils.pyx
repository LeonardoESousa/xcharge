import numpy as np
cimport cython

from libc.math cimport sqrt
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef distances(double[:] dx, double[:] dy, double[:] dz, int num):
    r = np.empty(num)
    cdef double [:] r_view = r
    cdef int i
    for i in range(num):
        r_view[i]  = sqrt(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i])     
    return r  

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef forster(double[:] Rf,int[:] mats,int num, double alpha_mu, double[:] r, double emi_rate):
    taxas = np.empty(num)
    cdef double [:] taxas_view = taxas 
    cdef double ratio
    cdef int i
    for i in range(num):
      if r[i] != 0:  
        ratio = Rf[mats[i]]/(alpha_mu+r[i])  
        ratio = ratio*ratio
        taxas_view[i] = ratio*ratio*ratio*emi_rate
      else:
        taxas_view[i] = 0.0
    return taxas  

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef forster_anni(double[:] Rf,int[:] mats,int num, double alpha_mu, double[:] r, double emi_rate, int[:] replace_pos, double[:] replace_raios, int mum):
    taxas = np.empty(num)
    cdef double [:] taxas_view = taxas 
    cdef double ratio
    cdef int [:] replace_view = replace_pos
    cdef int i
    for i in range(num):
      if r[i] != 0:
        ratio = Rf[mats[i]]/(alpha_mu+r[i])  
        ratio = ratio*ratio
        taxas_view[i] = ratio*ratio*ratio*emi_rate
      else:
        taxas_view[i] = 0.0
    for i in range(mum):
      ratio = replace_raios[i]/(alpha_mu+r[replace_view[i]])
      ratio = ratio*ratio
      taxas_view[replace_view[i]] = ratio*ratio*ratio*emi_rate   

    return taxas  

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef jump(double[:] jump_rate, int num, double random_number):
    cdef double soma
    cdef double [:] cumsum = np.empty(num)
    cdef double [:] jump_rate_view = jump_rate
    cdef int i
    for i in range(num):
        soma += jump_rate_view[i]
        cumsum[i] = soma
    random_number = random_number*soma
    for i in range(num):
        if random_number <= cumsum[i]:
            return soma,i

