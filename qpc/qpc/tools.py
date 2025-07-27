import numpy as np
from scipy.sparse import spdiags
from scipy.linalg import eigh
from scipy import sparse

import math
from numpy import linalg as LA
import scipy.sparse.linalg as sl

from scipy.sparse import coo_matrix

import random

from . _fastmul import fastmul
from . _fastmuli import fastmuli

from . _evolution import evolution

def mv(m, vin, vout, cin=1, cout=0):
    fastmul(m.data, m.indices, m.indptr, cin, vin, cout, vout)
    
def mvi(length, m, vin, vout, cin=1, cout=0):
    fastmuli(m.data, m.indices, m.indptr, cin, vin, cout, vout, length)

def tridiag(e, h):
    data = [np.concatenate((h, np.array([0]))), np.array(e), np.concatenate((np.array([0]), h))]
    diags = np.array([-1, 0, 1])

    n = len(e)
    return spdiags(data, diags, n, n).tocsc()

def dyad(ket, bra):
    return np.kron(ket, bra.T.conj())

def as_column_vector(ket):
    return ket[:, None]

def make_hermitean(m):
    m += m.T.conj()
    m /= 2

def find_largest_eigs(m, k = None):
    n = len(m)
    if k is None:
        k = n
    k = min(k, n)
    e, v = eigh(m, subset_by_index = [n - k, n - 1])
    e = np.flip(e)
    v = np.flip(v, axis = 1)
    return (e, v)

def find_smallest_eigs(m, k = None):
    n = len(m)
    if k is None:
        k = n
    k = min(k, n)
    e, v = eigh(m, subset_by_index = [0, k - 1])
    return (e, v)

def find_eigs_ascending(m):
    e, v = eigh(m)
    return (e, v)

def find_eigs_descending(m):
    e, v = eigh(m)
    e = np.flip(e)
    v = np.flip(v, axis = 1)
    return (e, v)
   
def evolutionpy_chained(dt, H, initial_state, start_time = None, start_index = None, end_time = None, end_index = None, tol = 10**(-6), first_in_chain = True):

    K = initial_state.size

    use_time = False
    use_index = False

    if not start_time is None and not end_time is None:

        use_time = True
        nt = math.floor((end_time - start_time) / dt)

    if not start_index is None and not end_index is None:

        use_index = True
        nt = end_index - start_index

    if not use_index != use_time:
        raise ValueError('evolution should be called either in time or in step-index mode')

    if (first_in_chain):

        if (use_index):
            yield (start_index, initial_state)
        else:
            yield (start_time, initial_state)

    psi = np.copy(initial_state)
    psi_mid = np.copy(psi)

    b = nt

    for i in range(0, b):

        psi_mid[:] = psi

        if (use_index):

            H_mid = H(start_index + i)

        else:

            H_mid = H(start_time + (i + 0.5) * dt)

        while(True):

            psi_mid_next = H_mid @ psi_mid

            psi_mid_next = psi - 1j * dt / 2 * psi_mid_next

            err = max(abs(psi_mid_next - psi_mid))

            swp = psi_mid_next
            psi_mid_next = psi_mid
            psi_mid = swp

            if err < tol:
                break

        psi = 2 * psi_mid - psi

        if (use_index):
            yield (start_index + i + 1, psi)
        else:
            yield (start_time + i * dt, psi)
    
def evolutionpy(dt, H, initial_state, start_time = None, start_index = None, end_time = None, end_index = None, tol = 10**(-6)):

    K = initial_state.size

    use_time = False
    use_index = False

    if not start_time is None and not end_time is None:

        use_time = True
        nt = math.floor((end_time - start_time) / dt)

    if not start_index is None and not end_index is None:

        use_index = True
        nt = end_index - start_index

    if not use_index != use_time:
        raise ValueError('evolution should be called either in time or in step-index mode')


    if (use_index):
        yield (start_index, initial_state)
    else:
        yield (start_time, initial_state)

    psi = np.copy(initial_state)
    psi_mid = np.copy(psi)

    b = nt - 1

    for i in range(0, b):

        psi_mid[:] = psi

        if (use_index):

            H_mid = H(start_index + i)

        else:

            H_mid = H(start_time + (i + 0.5) * dt)

        while(True):

            psi_mid_next = H_mid @ psi_mid

            psi_mid_next = psi - 1j * dt / 2 * psi_mid_next

            err = max(abs(psi_mid_next - psi_mid))

            swp = psi_mid_next
            psi_mid_next = psi_mid
            psi_mid = swp

            if err < tol:
                break

        psi = 2 * psi_mid - psi

        if (use_index):
            yield (start_index + i + 1, psi)
        else:
            yield (start_time + i * dt, psi)
            
def evolutionpy2(dt, apply_H, eval_O, initial_state, start_time = None, start_index = None, end_time = None, end_index = None, tol = 10**(-6)):

    K = initial_state.size

    use_time = False
    use_index = False

    if not start_time is None and not end_time is None:

        use_time = True
        nt = math.floor((end_time - start_time) / dt)

    if not start_index is None and not end_index is None:

        use_index = True
        nt = end_index - start_index

    if not use_index != use_time:
        raise ValueError('evolution should be called either in time or in step-index mode')


    if (use_index):
        eval_O(start_index, initial_state)
    else:
        eval_O(start_time, initial_state)

    psi = np.copy(initial_state)
    psi_mid = np.copy(psi)

    psi_mid_next = np.zeros(K, dtype = complex)

    for i in range(0, nt - 1):

        #psi_tmp = psi
        #psi = psi_mid
        #psi_mid = psi_tmp
        
        psi_mid[:] = psi

        while(True):

            psi_mid_next.fill(0)
            
            if (use_index):
                apply_H(start_index + i, psi_mid, psi_mid_next)
            else:
                apply_H(start_time + (i + 0.5) * dt, psi_mid, psi_mid_next)
            
            #psi_mid_next = H_mid @ psi_mid

            psi_mid_next = psi - 1j * dt / 2 * psi_mid_next

            err = max(abs(psi_mid_next - psi_mid))

            swp = psi_mid_next
            psi_mid_next = psi_mid
            psi_mid = swp

            if err < tol:
                break

        psi = 2 * psi_mid - psi

        if (use_index):
            eval_O(start_index + i + 1, psi)
        else:
            eval_O(start_time + i * dt, psi)

class LocalObservables:
    def __init__(self, f, m_max, n_max):
    
        #self.local_observables = self.local_projections(f, m_max, n_max)
        self.local_observables = self.local_projections_f(f, m_max, n_max)

    def local_projections_f(self, f, m_max, n_max):
        
        from qpc import secondquant
        from qpc import local_op
        from scipy.sparse import csc_matrix
        
        self.K = f.dimension
        self.local_dim = n_max + 1
        
        a = f.annihilate
        a_dag = f.create
        
        local_ops = []
        
        o_data = np.zeros(self.K, dtype = complex)
        o_ind = np.zeros(self.K, dtype = np.int32)
        
        o_ptr = np.zeros(self.K + 1, dtype = np.int32)
        for i in range(self.K + 1):
            o_ptr[i] = i
        
        for i in range(0, m_max + 1):
            
            o = np.array([f.occupations(j)[i] for j in range(self.K)])
            
            a_ = a[i]
            b_ = a_dag[i]
            
            mode_op = [[] for l in range(self.local_dim)]
            
            for p in range(self.local_dim):
                for q in range(self.local_dim):
                    
                    local_op(a_.data, a_.indices, a_.indptr, \
                        b_.data, b_.indices, b_.indptr, \
                            o_data, o_ind, o, p, q)
            
                    #print(o_data)
                    #print(o_ind)
                    #print(o_ptr)

                    mode_op[p].append(csc_matrix((np.copy(o_data), np.copy(o_ind), np.copy(o_ptr)), shape = (self.K, self.K)))
                    
            local_ops.append(mode_op)
            
        return local_ops
        
    def local_projections(self, f, m_max, n_max):
        # construct projections to local hilbertspaces

        self.K = f.dimension

        a = f.annihilate
        a_dag = f.create

        local_op = []

        self.local_dim = n_max + 1

        # for each mode

        for i in range(0, m_max + 1):
   
            mode_op = [ [ ([], [], []) for l in range(self.local_dim) ] for k in range(self.local_dim) ]  # ([data], [row], [col])

            # diagonal 

            for j in range(0, self.K):

                o = f.occupations(j)
                p = o[i]
                q = p

                mode_op[q][p][0].append(1.0)
                mode_op[q][p][1].append(j)
                mode_op[q][p][2].append(j)

            # upper diagonal

            c_dag = a_dag[i].tocoo()
            c_dag.eliminate_zeros()
            row = c_dag.row
            col = c_dag.col

            for j in col:

                o = f.occupations(j)
                p = o[i]
        
                j_ = j
                d_ = 1.0

                for q in range(p + 1, self.local_dim):

                    if (not j_ in col):
                        break

                    j_ = row[col.tolist().index(j_)]

                    o_ = f.occupations(j_)

                    assert q == o_[i]


                    mode_op[q][p][0].append(d_)
                    mode_op[q][p][1].append(j_)
                    mode_op[q][p][2].append(j)

            # lower diagonal

            c = a[i].tocoo()
            c.eliminate_zeros()
            row = c.row
            col = c.col

            for j in col:

                o = f.occupations(j)
                p = o[i]
        
                j_ = j
                d_ = 1.0

                for q in reversed(range(0, p)):

                    if (not j_ in col):
                        break

                    j_ = row[col.tolist().index(j_)]

                    o_ = f.occupations(j_)
                    q = o_[i]

                    mode_op[q][p][0].append(d_)
                    mode_op[q][p][1].append(j_)
                    mode_op[q][p][2].append(j)

            for k in range(self.local_dim):
                for l in range(self.local_dim):
                    s = mode_op[k][l]
                    mode_op[k][l] = coo_matrix((s[0], (s[1], s[2])), shape = (self.K, self.K), dtype = complex).tocsc()

            local_op.append(mode_op)

        return local_op

    def partial_trace(self, psi, measured_mode):

        rho = np.zeros((self.local_dim, self.local_dim), dtype = complex)

        for i in range(self.local_dim):
            for j in range(self.local_dim):
                rho[i, j] = np.vdot(psi, self.local_observables[measured_mode][j][i] @ psi)

        return rho

    def project_to_vacuum(self, psi, measured_mode):

        psi_vac = self.local_observables[measured_mode][0][0] @ psi  
        return psi_vac

    def quantum_jump(self, psi, measured_mode):

        # find the "preferred" jump basis
        rho = self.partial_trace(psi, measured_mode)
        jump_probs, jump_states = LA.eigh(rho)

        # select the jump 

        xi = random.random()

        jump_index = 0

        for k in range(0, self.local_dim):
            xi = xi - jump_probs[k]

            if xi < 0:
                jump_index = k
                break

        # execute the jump

        psi_collapsed = np.zeros(self.K, dtype = complex)

        for k in range(self.local_dim):
            psi_collapsed = psi_collapsed + np.conj(jump_states[k, jump_index]) * self.local_observables[measured_mode][0][k] @ psi  

        # normalize
            
        psi_collapsed = psi_collapsed / np.sqrt(np.vdot(psi_collapsed, psi_collapsed))

        return psi_collapsed
    
    def quantum_jumpEx(self, psi, measured_mode):

        # find the "preferred" jump basis
        rho = self.partial_trace(psi, measured_mode)
        jump_probs, jump_states = LA.eigh(rho)

        # select the jump 

        xi = random.random()

        jump_index = 0

        for k in range(0, self.local_dim):
            xi = xi - jump_probs[k]

            if xi < 0:
                jump_index = k
                break

        # execute the jump

        psi_collapsed = np.zeros(self.K, dtype = complex)

        for k in range(self.local_dim):
            psi_collapsed = psi_collapsed + np.conj(jump_states[k, jump_index]) * self.local_observables[measured_mode][0][k] @ psi  

        # normalize
            
        psi_collapsed = psi_collapsed / np.sqrt(np.vdot(psi_collapsed, psi_collapsed))

        return (psi_collapsed, jump_probs, jump_states, jump_index)
    
    def a(self):
        
        a_ = np.zeros((self.local_dim, self.local_dim), dtype = complex)
        
        for i in range(self.local_dim-1):
            a_[i, i + 1] = np.sqrt(i+1)
          
        return a_
        
    def a_dag(self):
        
        a_dag_ = np.zeros((self.local_dim, self.local_dim), dtype = complex)
        
        for i in range(self.local_dim-1):
            a_dag_[i + 1, i] = np.sqrt(i+1)
          
        return a_dag_

    def pair_quantum_jump(self, psi_r, psi_l, measured_mode):

        rho_l = self.partial_trace(psi_l, measured_mode)
        c_l, phi_l = LA.eigh(rho_l)

        c_r = np.zeros(self.local_dim, dtype = np.cdouble)

        psi_r_collapsed = np.zeros((self.K, self.local_dim), dtype = complex)

        for k in range(self.local_dim):

            for l in range(self.local_dim):
                psi_r_collapsed[:, k] = psi_r_collapsed[:, k] + np.conj(phi_l[l, k]) * self.local_observables[measured_mode][0][l] @ psi_r  

            norm = np.vdot(psi_r_collapsed[:, k], psi_r_collapsed[:, k])

            if norm > 10**(-9):
                c_r[k] = math.sqrt(norm)
                psi_r_collapsed[:, k] = psi_r_collapsed[:, k] / c_r[k]
            else:
                c_r[k] = 0

        jump_probs = c_l * c_r
        z = np.sum(jump_probs)
        jump_probs = jump_probs / z

        # select the jump 
        xi = random.random()
        jump_index = 0

        for k in range(0, self.local_dim):
            xi = xi - jump_probs[k]

            if xi < 0:
                jump_index = k
                break

        psi_l_collapsed = np.zeros(self.K, dtype = complex)

        for k in range(self.local_dim):
            psi_l_collapsed = psi_l_collapsed + np.conj(phi_l[k, jump_index]) * self.local_observables[measured_mode][0][k] @ psi_l  

        # normalize
            
        psi_l_collapsed = psi_l_collapsed / np.sqrt(np.vdot(psi_l_collapsed, psi_l_collapsed))

        return (psi_r_collapsed[:, jump_index], psi_l_collapsed, z)

class spin_boson_model:
    def __init__(self, num_modes, max_num_quanta):
        
        from qpc import secondquant as sq
    
        hs_atom = sq.fock_space(num_modes = 1, max_total_occupation = 1, statistics = 'Bose')
        m = num_modes
        n = max_num_quanta
        fs_chain = sq.fock_space(num_modes = m, max_total_occupation = n, statistics = 'Bose') 
        hs_joint = sq.fock_space_kron(hs_atom, fs_chain)
        b_hat = hs_joint.annihilate
        b_hat_dag = hs_joint.create
        sigma_m = b_hat[0]
        sigma_p = b_hat_dag[0]
        sigma = [hs_joint.sigmax(0), hs_joint.sigmay(0), hs_joint.sigmaz(0)]
        a_hat = b_hat[1:]
        a_hat_dag = b_hat_dag[1:]
    
        self.hs_joint = hs_joint
        
        self.s_m = sigma_m
        self.s_p = sigma_p
        self.s = sigma
        
        self.s_x = sigma[0]
        self.s_y = sigma[1]
        self.s_z = sigma[2]
                         
        self.a = a_hat
        self.a_dag = a_hat_dag
        
        self.eye = hs_joint.eye
        
        self.num_modes = num_modes
        self.max_num_quanta = max_num_quanta
        
        self.dimension = hs_joint.dimension
        
    def get_local_observables(self):
        return LocalObservables(self.hs_joint, self.num_modes, self.max_num_quanta)
    
class fermion_2lead_fermion_model:
    def __init__(self, num_impurity_modes, num_reservoir_modes, max_num_quanta):
        
        from qpc import secondquant as sq
        
        self.n_qua = max_num_quanta
        self.m_imp = num_impurity_modes
        self.m_env = num_reservoir_modes
         
        
        total_num_modes = num_impurity_modes + 2 * num_reservoir_modes
        
        self.m_tot = total_num_modes
        
        joint = sq.fock_space(num_modes = total_num_modes, max_total_occupation = max_num_quanta, statistics = 'Fermi')
        
        self.d = joint.annihilate[: num_impurity_modes] # impurity modes
        self.d_dag = joint.create[: num_impurity_modes]
        
        self.l = joint.annihilate[num_impurity_modes : num_impurity_modes + num_reservoir_modes] # reservoir modes
        self.l_dag = joint.create[num_impurity_modes : num_impurity_modes + num_reservoir_modes]
        
        self.r = joint.annihilate[num_impurity_modes + num_reservoir_modes :] # reservoir modes
        self.r_dag = joint.create[num_impurity_modes + num_reservoir_modes :]
        
        self.space = joint
        self.dimension = joint.dimension
        
    def get_local_observables(self):
        return LocalObservables(self.space, self.m_tot - 1, 1)
    
class fermion_fermion_model:
    def __init__(self, num_impurity_modes, num_reservoir_modes, max_num_quanta):
        
        from qpc import secondquant as sq
        
        self.n_qua = max_num_quanta
        self.m_imp = num_impurity_modes
        self.m_env = num_reservoir_modes
         
        
        total_num_modes = num_impurity_modes + num_reservoir_modes
        
        self.m_tot = total_num_modes
        
        joint = sq.fock_space(num_modes = total_num_modes, max_total_occupation = max_num_quanta, statistics = 'Fermi')
        
        self.d = joint.annihilate[: num_impurity_modes] # impurity modes
        self.d_dag = joint.create[: num_impurity_modes]
        
        self.c = joint.annihilate[num_impurity_modes :] # reservoir modes
        self.c_dag = joint.create[num_impurity_modes :]
        
        self.space = joint
        self.dimension = joint.dimension
        
    def get_local_observables(self):
        return LocalObservables(self.space, self.m_tot - 1, 1)
        
class top_boson_model:
    def __init__(self, j, num_modes, max_num_quanta):
        
        from qpc import secondquant as sq
    
        top = sq.fock_space(num_modes = 1, max_total_occupation = round(2*j), statistics = 'Bose')
    
       # hs_atom = sq.fock_space(num_modes = 1, max_total_occupation = 1, statistics = 'Bose')
        m = num_modes
        n = max_num_quanta
        fs_chain = sq.fock_space(num_modes = m, max_total_occupation = n, statistics = 'Bose') 
        hs_joint = sq.fock_space_kron(top, fs_chain)
        b_hat = hs_joint.annihilate
        b_hat_dag = hs_joint.create
        #sigma_m = b_hat[0]
        #sigma_p = b_hat_dag[0]
        
        j_m = hs_joint.j_m[0]
        j_p = hs_joint.j_p[0]
        j_x = hs_joint.j_x[0]
        j_y = hs_joint.j_y[0]
        j_z = hs_joint.j_z[0]
        
        #sigma = [hs_joint.sigmax(0), hs_joint.sigmay(0), hs_joint.sigmaz(0)]
        a_hat = b_hat[1:]
        a_hat_dag = b_hat_dag[1:]
    
        self.hs_joint = hs_joint
        
        self.j_m = j_m
        self.j_p = j_p
        self.j_x = j_x
        self.j_y = j_y
        self.j_z = j_z
        self.j = j
                         
        self.a = a_hat
        self.a_dag = a_hat_dag
        
        self.eye = hs_joint.eye
        
        self.num_modes = num_modes
        self.max_num_quanta = max_num_quanta
        
        self.dimension = hs_joint.dimension
        
    def kick_z(self, a, b):
        
        j_z = self.hs_joint.f1.j_z[0]
        
        eye_1 = self.hs_joint.f1.eye
        eye_2 = self.hs_joint.f2.eye
        
        #sparse.kron(f1.create[k], f2.eye).tocsc()
        
        return sparse.kron(sl.expm(-1j * a * (j_z - b * eye_1) @ (j_z - b * eye_1)), eye_2).tocsc()
        
    def get_local_observables(self):
        return LocalObservables(self.hs_joint, self.num_modes, self.max_num_quanta)
    
def forward_spatial_lightcone_generator(chain_generator, dt, chunk_size=100, guard_size=10, time_chunk = None, rcut = 10**(-6)):
        
        
    e = []
    h = []
    
    if chunk_size < 2:
        raise Exception('chunk_size must be >= 2')

    if time_chunk is None:
        time_chunk = 100 * dt
        
    for i in range(chunk_size):
        try:
            e_, h_ = next(chain_generator)
        except StopIteration:
            break
            
        e.append(e_)
        h.append(h_)
  
    ns = len(e)

    H1 = tridiag(e, h[:-1])
        
    phi_ini = np.zeros(ns, dtype = np.cdouble)
    phi_ini[0] = 1
        
    m = 1
    a = 0
        
    rho_plus = np.zeros((ns, ns), dtype = complex)
        
    n_time_chunk = round(time_chunk / dt)
        
    a_chunk = 0
    b_chunk = n_time_chunk
    phi_begin_chunk = phi_ini
    phi_begin = np.copy(phi_ini)
        
    first_chunk = True
    first_interval = True

    while(True):
                
        end_of_chain = False
        i_end = b_chunk
        phi_end = None
        
    
        def H1t(ti):
            return H1
        
        for i, phi in evolutionpy_chained(start_index = a_chunk, end_index = b_chunk, H = H1t, dt = dt, initial_state = phi_begin_chunk, first_in_chain = first_chunk):
                
            if not end_of_chain:
                
                phi_ = as_column_vector(phi)
                rho_plus += dyad(phi_, phi_) * dt
                    
                # find max eigenvalue
                pi_max, _ = find_largest_eigs(rho_plus, 1)
                    
                # check whether the next site is inside the lightcone
                site_sig = rho_plus[m + 1 - 1, m + 1 - 1]
                lr_metric = site_sig - rcut * pi_max
                if lr_metric > 0:
                        
                    b = i
                    #yield (a, b, m, evolutionpy_chained(start_index = a, end_index = b, H = H1t, dt = dt, initial_state = phi_begin, first_in_chain = first_interval))
                    yield (a, b, m, evolutionpy_chained(start_index = a, end_index = b - 1, H = H1t, dt = dt, initial_state = phi_begin, first_in_chain = True))
                         
                    first_interval = False

                    m += 1
                    a = i
                    phi_begin = np.copy(phi)
                        
                    # need to enlong the chain
                    if m > ns - guard_size:
                            
                        end_of_chain = True
                        i_end = i
                        phi_end = np.copy(phi)

        first_chunk = False
                            
        if end_of_chain:
                                
            for i in range(chunk_size):
                try:
                    e_, h_ = next(chain_generator)
                except StopIteration:
                    break
                    
                e.append(e_)
                h.append(h_)
                    
            H1 = tridiag(e, h[:-1])
            
            ns_ = len(e)
            
            phi_begin = np.concatenate((phi_begin, np.zeros(ns_ - ns, dtype = complex)))
            
            rho_plus_ = np.zeros((ns_, ns_), dtype = complex)
            rho_plus_[0:ns, 0:ns] = rho_plus
                
            
                
            b_chunk = i_end
            phi = phi_end
                
            phi = np.concatenate((phi, np.zeros(ns_ - ns, dtype = complex)))

            ns = ns_
            rho_plus = rho_plus_

        a_chunk = b_chunk
        b_chunk += n_time_chunk
            
        phi_begin_chunk = phi
                
        
        




