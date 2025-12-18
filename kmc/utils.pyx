import numpy as np
cimport cython
cimport numpy as cnp
cnp.import_array()

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
    cdef double soma = 0.0
    if num == 0:
        return 0.0, -1
    cdef double[:] cumsum = np.empty(num)
    cdef double[:] jump_rate_view = jump_rate
    cdef int i
    for i in range(num):
        soma += jump_rate_view[i]
        cumsum[i] = soma
    if soma <= 0:
        return 0.0, -1
    random_number = random_number * soma
    for i in range(num):
        if random_number <= cumsum[i]:
            return soma, i


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef tuple select_event(
    double[:] dx,
    double[:] dy,
    double[:] dz,
    int[:] mats,
    int matlocal,
    object hops,
    object monos,
    object particle,
    object system,
    object cut,
    double randu):
    """
    Hot-path helper used by kmc.__main__.decision.
    Returns (total_rate, chosen_index, sizes_list).
    """
    # Ensure downstream Python rate functions see NumPy arrays (not memoryviews)
    cdef cnp.ndarray[cnp.double_t, ndim=1] dx_arr = np.asarray(dx)
    cdef cnp.ndarray[cnp.double_t, ndim=1] dy_arr = np.asarray(dy)
    cdef cnp.ndarray[cnp.double_t, ndim=1] dz_arr = np.asarray(dz)
    cdef cnp.ndarray[cnp.int_t, ndim=1] mats_arr = np.asarray(mats)

    cdef Py_ssize_t n = dx_arr.shape[0]
    cdef cnp.ndarray[cnp.double_t, ndim=1] r = distances(dx_arr, dy_arr, dz_arr, n)

    cdef list rates = []
    cdef list sizes = []
    cdef object proc
    cdef object arr_any

    for proc in hops:
        arr_any = proc.rate(r=r, dx=dx_arr, dy=dy_arr, dz=dz_arr, system=system, particle=particle, mats=mats_arr, matlocal=matlocal, cut=cut)
        rates.append(arr_any)
        sizes.append(len(arr_any))

    for proc in monos:
        arr_any = [proc.rate(material=matlocal)]
        rates.append(arr_any)
        sizes.append(1)

    cdef Py_ssize_t total_size = 0
    cdef Py_ssize_t size_val
    for size_val in sizes:
        total_size += size_val

    if total_size == 0:
        return (0.0, -1, sizes)

    cdef cnp.ndarray[cnp.double_t, ndim=1] flat = np.empty(total_size, dtype=np.double)
    cdef double[:] flat_view = flat
    cdef Py_ssize_t offset = 0
    cdef Py_ssize_t i
    for arr_any in rates:
        for i in range(len(arr_any)):
            flat_view[offset + i] = arr_any[i]
        offset += len(arr_any)

    cdef double soma
    cdef int idx
    soma, idx = jump(flat_view, total_size, randu)
    return (soma, idx, sizes)


from libc.math cimport exp
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef dexter(double[:] Rd,double invL,double emi_rate,int[:] mats,int num, double[:] r):
    taxas = np.empty(num)
    cdef double [:] taxas_view = taxas 
    cdef double ratio
    cdef int i
    for i in range(num):
      if r[i] != 0:  
        ratio = Rd[mats[i]]*invL
        taxas_view[i] = ratio*ratio*exp(2*invL*(Rd[mats[i]]-r[i]))*emi_rate
      else:
        taxas_view[i] = 0.0
    return taxas
