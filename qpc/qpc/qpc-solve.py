#!/usr/bin/env python3

import os
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import sys
import numpy as np
import qpc.tools as tools
from qpc.tools import mv
import time
from qpc import evolution_chained2_kicked
import math
import random
import pathlib



print("Wellcome to lightcone impurity solver version 10.03.2023-11.51")
print("Will read the lightcone info from the current directory")

#### execution parameters

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--fout', type=str, default = "results")
parser.add_argument('--fin', type=str, default = "lc")
parser.add_argument('--cfrom', type=int, default = 1)
parser.add_argument('--cto', type=int, default = 10)
parser.add_argument('--csize', type=int, default = 10)
parser.add_argument('--maxcores', type=int, default = 1)
parser.add_argument('--nquanta', type=int, default = 3)
parser.add_argument('--e', type=int, default = 1)
parser.add_argument('--h', type=int, default = 0.05)
parser.add_argument('--e_d', type=float, default = 1.0)
parser.add_argument('--U', type=int, default = 0.0)
parser.add_argument('--tmax', type=float, default = 200)


args = parser.parse_args()

chunk_from = args.cfrom
chunk_to = args.cto # inclusive
trajectories_per_chunk = args.csize
max_cores_used = args.maxcores

print("Will compute ", (chunk_to - chunk_from + 1) * trajectories_per_chunk, " trajectories")
print("on max ", max_cores_used, " CPU cores")

#### where to save the results

folder = args.fout

print("Will save results into ", folder)

#### where to find lightcones data

fin = args.fin
fin = os.path.join(os.getcwd(), fin)

print("Will load lightcones from ", fin)

#### Max simulation time and time step

tmax = args.tmax

print("Will simulate up to time", tmax)

####

    
#### Max number of quanta

n_max = args.nquanta

print("Max number of coupled quanta: ", n_max)


####

e_chain = args.e

print("on-site energy: ", e_chain)

####

h_chain = args.h

print("hopping: ", h_chain)
                    
####

e_d = args.e_d

print("Gate voltage e_d: ", e_d)
                    
####

U = args.U

print("Coulomb interaction U: ", U)
                    
####

#### Then the code follows

#### Import the lightcone info

#----------------------------------------------------

intervals = []
p = os.path.join(fin, "intervals.txt")
try:
    with pathlib.Path(p).open() as f:
        while True:
            l = [int(e) for e in next(f).split()] 
            intervals.append(l)
except StopIteration as e:
    pass

#----------------------------------------------------

couplings = []
p = os.path.join(fin, "couplings.txt")
try:
    with pathlib.Path(p).open() as f:
        while True:
            re = np.asarray([float(e) for e in next(f).split()])
            im = np.asarray([float(e) for e in next(f).split()])
            couplings.append(re + 1j * im)
except StopIteration as e:
    pass

#----------------------------------------------------

rotations = []

m_max = 0
n_rel = 0

p = os.path.join(fin, "rotations.txt")
with pathlib.Path(p).open() as f:

    for i in intervals:

        a = i[2]
        b = i[3]

        n_rel = max(n_rel, b)
        m_max = max(m_max, b - a)

        re = np.zeros((b - a, b - a), dtype = np.cdouble)
        for p in range(b - a):
            re[p, :] = np.asarray([float(e) for e in next(f).split()])

        im = np.zeros((b - a, b - a), dtype = np.cdouble)
        for p in range(b - a):
            im[p, :] = np.asarray([float(e) for e in next(f).split()])

        w = re + 1j * im

        rotations.append(w)

print("Max number of coupled modes: ", m_max)

#----------------------------------------------------

p = os.path.join(fin, "time.txt")
with pathlib.Path(p).open() as f:
    tg = np.loadtxt(f)
ntg = tg.size
dt = tg[1] - tg[0]

#-----------------------------------------------------

t = tg
nt = ntg
nt_ = np.where(t <= tmax)[0][-1]
t_ = t[0 : nt_]

    
#-----------------------------------------------------

p = os.path.join(os.getcwd(), folder, "time.txt")
os.makedirs(os.path.dirname(p), exist_ok = True)

with open(p, "w") as f:
    for i in range(nt_):
        print(str(t[i]), file = f)

#-----------------------------------------------------

print('computing sparse matrices...', flush = True)

m = tools.fermion_2lead_fermion_model(num_impurity_modes = 2, num_reservoir_modes = m_max, max_num_quanta = n_max)

lo = m.get_local_observables()

#-----------------------------------------------------

print('...done', flush = True)

#### Initial condition

psi_ini = np.zeros(m.space.dimension, dtype = complex)
psi_ini[0] = 1                    
                    
####

to_ring = [ _ % m_max for _ in range(0, n_rel)]

###

def Hdot(ti):

    # gate voltage 
    Hgate = e_d * m.d_dag[1] @ m.d[1] - e_d * m.d_dag[0] @ m.d[0]

    # coulomb interaction
    Hcoulomb = U * m.d_dag[1] @ m.d[1] - U * m.d_dag[0] @ m.d[0] @ m.d_dag[1] @ m.d[1]

    # hopping between the two dots
    Hhopping = h_chain * m.d_dag[1] @ m.d_dag[0] + h_chain *  m.d[0] @ m.d[1]

    return Hgate + Hcoulomb + Hhopping

la = lo.a()
la_dag = lo.a_dag()

####

def job(image):
    
    print('running chunk ', image, ', of ', trajectories_per_chunk, ' trajectories', flush = True)
    
    start_time = time.time()
    
    seed = 1000 * math.pi * image + pow(image, 2) 
    random.seed(seed)
    
    nqx = np.zeros(nt_, dtype = complex)
    #obs = sum([m.d_dag[i] @ m.d[i] for i in range(m.m_imp)])
    
    #j_probs = np.zeros((nt_, lo.local_dim))
    #j_disps = np.zeros((nt_, lo.local_dim))
    #j_ns    = np.zeros((nt_, lo.local_dim))

    j_probs = []
    j_disps = []
    j_ns = []
    
    j_av_disp = []
    
    psi = np.zeros(m.dimension, dtype = complex)
    psi_mid = np.zeros(m.dimension, dtype = complex)
    psi_mid_next = np.zeros(m.dimension, dtype = complex)
    psi_buff = np.zeros(m.dimension, dtype = complex)
        
    ####
    
    job.n_out = 0j
    job.obs = None

    def eval_o(ti, psi):
       
        mv(job.obs, psi, psi_buff)
        nqx[ti] = nqx[ti] + job.n_out + np.vdot(psi, psi_buff).real
        
    ####
    
    def eval():
        
        job.n_out = 0j
        job.obs = None
        
        psi_begin = np.copy(psi_ini)
    
        first_in_chain = 1

        ni = 0

        for (i, i_) in zip(intervals, [None] + intervals[:-1]):

            if i[0] >= nt_:
                return
            
            i1 = min(i[1], nt_-1)

            b = i[3]
            a = i[2]

            if not i_ is None:
                a_ = i_[2]
                for q in range(a_, a):
                    measured_l_mode = to_ring[q] + m.m_imp
                    psi_begin, j_l_prob, j_l_psi, j_ind  = lo.quantum_jumpEx(psi_begin, measured_l_mode)
                    
                    measured_r_mode = to_ring[q] + m.m_imp + m.m_env
                    psi_begin, j_r_prob, j_r_psi, j_ind  = lo.quantum_jumpEx(psi_begin, measured_r_mode)
                    
                    #j_probs.append(j_prob)
                    
                    #dd = np.zeros(j_psi.shape[1], dtype = complex)
                    
                    #for r in range(j_psi.shape[1]):
                    #    dd[r] = np.vdot(j_psi[:, r], la @ j_psi[:, r])
                        
                    #j_disps.append(dd)
                    
                    #ad = 0
                    
                    #for r in range(j_l_psi.shape[1]):
                    #    ad = ad + j_l_prob[r] * np.vdot(j_l_psi[:, r], la_dag @ la @ j_l_psi[:, r])
                    
                    #job.n_out = job.n_out + ad
                    
                    ad = np.vdot(j_r_psi[:, j_ind], la_dag @ la @ j_r_psi[:, j_ind])
                    
                    job.n_out = job.n_out + ad
                    
                    #j_av_disp.append(ad)
                    
            w = rotations[ni]

            l_ring = [m.l[to_ring[_]] for _ in range(a, b)]
            l_ring_dag = [m.l_dag[to_ring[_]] for _ in range(a, b)]
            
            r_ring = [m.r[to_ring[_]] for _ in range(a, b)]
            r_ring_dag = [m.r_dag[to_ring[_]] for _ in range(a, b)]

            Hw = m.space.emptyH
            if not w is None:
                for p in range(b - a):
                    for q in range(b - a):
                        Hw += 1j * r_ring_dag[q] @ r_ring[p] * w[q, p].conj()
                        Hw += 1j * l_ring_dag[q] @ l_ring[p] * w[q, p]
                        
            job.obs = m.space.emptyH
            for p in range(b - a):
                job.obs = job.obs + r_ring_dag[p] @ r_ring[p]

            def Vint(ti):
                V_r = m.d_dag[1] @ sum(couplings[ti] * r_ring)
                V_l = - sum(couplings[ti] * l_ring_dag) @ m.d[0]
                V = V_r + V_l
                V = V + V.conj().transpose()
                return(V)            
          
            def Ht(ti):
                return Hdot(ti) + Vint(ti) + Hw

            def Hwt(ti):
                return Hw

            eval.Ht_ = None

            def begin_step(ti, psi):
                eval.Ht_ = Ht(ti)

            def apply_H(ti, psi_in, psi_out):
                mv(eval.Ht_, psi_in, psi_out)

            evolution_chained2_kicked(i[0], i1, dt, begin_step, apply_H, eval_o, psi_begin, psi, psi_mid, psi_mid_next, first_in_chain)
            first_in_chain = 0

            psi_begin = np.copy(psi)

            ni = ni + 1
           
    for i in range(trajectories_per_chunk):
        eval()

    end_time = time.time()
    
    print('chunk ', image, ' execution time: ', end_time-start_time, ' sec', flush = True)

    ####
    
    path = os.path.join(os.getcwd(), folder, "x_" + str(image) + ".txt")
    
    with open(path, "w") as f:
        print(str(trajectories_per_chunk) + " " + str(0), file = f)
        for i in range(nt_):
            print(str(nqx[i].real) + " " + str(nqx[i].imag), file = f)
                    

def job_(image):
    try:
        job(image)
    except:
        import traceback
        traceback.print_exception(*sys.exc_info())
        raise()

def progress_indicator(f):
    global progress
    progress += 1
    print('finished ', progress, 'job  out of ', chunk_to - chunk_from + 1, flush = True)

if __name__ == '__main__':
    
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures import wait
        
    with ProcessPoolExecutor(max_workers = max_cores_used) as executor:
        
        print('submitting jobs to pool ...', flush = True)
        
        futures = [executor.submit(job_, i) for i in range(chunk_from, chunk_to + 1)]
        global progress
        progress = 0
        for f in futures:
            f.add_done_callback(progress_indicator)
            
        wait(futures)
        
        for f in futures:
            exception = f.exception()
            if not exception is None:
                import traceback
                print("Exception:")
                print(exception)
                traceback.print_exc()
            
    print('\nDone!')
