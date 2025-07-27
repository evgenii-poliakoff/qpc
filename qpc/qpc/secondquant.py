    #!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy import linalg as LA
from scipy import sparse
from scipy.sparse import coo_matrix
from scipy.sparse import coo_array
from scipy.sparse import identity
from scipy.sparse import diags



class fock_space:
    def __init__(self, num_modes, max_total_occupation, statistics, fermi_sea_modes = None, max_local_occupations = None):

        self.global_exc = max_total_occupation #internal parameter
        self.statistics = statistics #internal parameter
            
        if statistics == 'Bose':
            
            if max_local_occupations is None:
        #1)
                self.modes = num_modes
        #2)
                self.local_exc = np.full(self.modes, max_total_occupation)
            else:
        #1)
                self.modes = num_modes #len(about_excitations)
        #2)
                self.local_exc = np.array(max_local_occupations)
            self.local_exc1 = np.array(self.local_exc+1) #internal parameter
        
        elif statistics == 'Fermi':
            if max_local_occupations is None:
        #1)
                self.modes = num_modes
        #2)
                
            else:
        #1)
                self.modes = num_modes
        #2)
                
        
        
        elif statistics == 'Fermi_sea':
            if max_local_occupations is None:
        #1)
                self.modes = num_modes
        #2)
                
            else:
        #1)
                self.modes = len(about_excitations)
        #2)
                
            
            self.fermi_sea_modes = fermi_sea_modes
                
                
        else:
            print('РЎhoose the statistics: Bose, Fermi or Fermi_sea') 
            
                    
        if not max_local_occupations is None and not len(max_local_occupations) == num_modes:
            raise Exception("Number of modes should be equal to number of local occupation constrains")
        
        
        
        #3)
        if statistics == 'Fermi_sea':
            self.states_list = list(self.sea_generator())
        else:
            self.states_list = list(self.states_generator())
        
        #4)
        self.find_index = {state: idx for idx, state in enumerate(self.states_list)}
        
        #5)
        self.dimension  = len(self.states_list)
        
        #6)
        self.emptyH = coo_matrix((self.dimension , self.dimension ), dtype = complex).tocsc()
        
        #7)
        self.eye = sparse.eye(self.dimension ).tocsc()
        
        
        if statistics == 'Bose':
            
        #8)    
            self.annihilate=[]
            current_state = []
            for k in range(self.modes):
                row = np.zeros(self.dimension )
                col = np.arange(self.dimension , dtype=int)
                data = np.zeros(self.dimension )
                for i in range(self.dimension ):
                    if self.states_list[i][k]==0:
                        row[i] = i
                    else:
                        current_state = list(self.states_list[i])
                        current_state[k]-=1
                        p = self.index(current_state)
                        row[i] = p
                        data[i]=(self.states_list[i][k])**0.5
                A = coo_matrix((data, (row, col)), shape=(self.dimension , self.dimension ), dtype = complex)
                #A.eliminate_zeros()
                self.annihilate.append(A.tocsc()) 
                
            
        #9)
            self.create=[]
            current_state = []
            for k in range(self.modes):
                row = np.zeros(self.dimension )
                col = np.arange(self.dimension , dtype=int)
                data = np.zeros(self.dimension )
                for i in range(self.dimension ):
                    if self.states_list[i][k] == self.local_exc[k] or sum(self.states_list[i]) == self.global_exc:#!!!!!!!!!!!!!!!!!!!!!!
                        row[i] = i 
                    else:
                        current_state = list(self.states_list[i])
                        current_state[k]+=1
                        p = self.index(current_state)
                        row[i] = p
                        data[i]=(self.states_list[i][k]+1)**0.5
                A = coo_matrix((data, (row, col)), shape=(self.dimension , self.dimension ), dtype = complex)
                #A.eliminate_zeros()
                self.create.append(A.tocsc())
       
        elif statistics == 'Fermi' or 'Fermi_sea':
        
        #8)
            self.annihilate=[]
            current_state = []
            for k in range(self.modes):
                row = np.zeros(self.dimension )
                col = np.arange(self.dimension , dtype=int)
                data = np.zeros(self.dimension )
                for i in range(self.dimension ):
                    if self.states_list[i][k]==0:
                        row[i] = i
                    else:
                        current_state = list(self.states_list[i])
                        y = sum(current_state[:k])
                        current_state[k] = 0
                        p = self.index(current_state)
                        
                        if (p=='666'):
                            continue
                        row[i] = p
                        data[i]=(-1)**y*(self.states_list[i][k])**0.5
                A = coo_matrix((data, (row, col)), shape=(self.dimension , self.dimension ), dtype = complex)
                #A.eliminate_zeros()
                self.annihilate.append(A.tocsc())
                              
         #9)       
            self.create=[]
            current_state = []
            for k in range(self.modes):
                row = np.zeros(self.dimension )
                col = np.arange(self.dimension , dtype=int)
                data = np.zeros(self.dimension )
                for i in range(self.dimension ):
                    if self.states_list[i][k] == 1:
                        row[i] = i 
                    else:
                        current_state = list(self.states_list[i])
                        y = sum(current_state[:k])
                        current_state[k] = 1
                        p = self.index(current_state)
                        if (p=='666'):
                            continue
                        row[i] = p
                        data[i]=(-1)**y*(self.states_list[i][k]+1)**0.5
                A = coo_matrix((data, (row, col)), shape=(self.dimension , self.dimension ), dtype = complex)
                #A.eliminate_zeros()
                self.create.append(A.tocsc())
            
        else:
            print('Сhoose the statistics: Bose, Fermi or Fermi_sea')
                  
                
    def sigmax(self, i):

        return(self.annihilate[i] + self.create[i])
    
    
    def sigmay(self, i):
 
        return(-1j*(-self.annihilate[i] + self.create[i]))
    
    
    def sigmaz(self, i):

        return(- self.annihilate[i]@self.create[i] + self.create[i]@self.annihilate[i])
    
    
    def occupations(self,i):
     
        if i >= self.dimension :
            print('the number is out of range')
            
        else:
            return(np.array(self.states_list[i]))
        
    
    
    def index(self, state):
       
        if len(state) != self.modes:
            print('incorrect size of an array')
        else:
            s = tuple(state)
            try:
                return(self.find_index[s])
            except:
                return('666')
    
    def states_generator(self):
      
        if self.statistics == 'Bose':
            #Generates all the possible states for given Fock space
            current_state = (0,)*len(self.local_exc1)
            n = 0
            while True:
                yield current_state
                j = len(self.local_exc1) - 1
                current_state = current_state[:j] + (current_state[j]+1,)
                n += 1
                while n > self.global_exc or current_state[j] >= self.local_exc1[j]:
                    j -= 1
                    if j < 0:
                        return
                    n -= current_state[j+1] - 1
                    current_state = current_state[:j] + (current_state[j]+1, 0) + current_state[j+2:]
      
        elif self.statistics == 'Fermi':
            modes = self.modes
            #Generates all the possible states for given Fock space
            current_state = (0,)*modes
            n = 0
            while True:
                yield current_state
                j = modes - 1
                current_state = current_state[:j] + (current_state[j]+1,)
                n += 1
                while n > self.global_exc or current_state[j] >= 2:
                    j -= 1
                    if j < 0:
                        return
                    n -= current_state[j+1] - 1
                    current_state = current_state[:j] + (current_state[j]+1, 0) + current_state[j+2:]          
        
    def occupations(self,i):
     
        if i >= self.dimension :
            print('the number is out of range')
            
        else:
            return(np.array(self.states_list[i]))
      

    def sea_generator(self):
        m_l = self.modes
        mp_l = self.fermi_sea_modes
        mh_l = m_l - mp_l
        N = self.global_exc
        if N>(m_l):
            N = m_l
        for Nl in range(N+1):
            Fl1 = min(Nl,mp_l)
            Fl2 = min(Nl,mh_l)                
            for i in range (Nl-Fl1,Fl2+1):
                for hl1 in self.generator_particles(mp_l, Nl-i):    
                    for hl2 in self.generator_holes(mh_l,i):
                        yield (hl1 + hl2)

    def generator_particles(self,m,N):
        if (m == 0):
            yield ()
            return
        current_state = (1,)*m
        n = 0
        while True:
            if (n == N):
                yield current_state
            j = 0
            current_state = (current_state[j]-1,) + current_state[j+1:]
            n += 1
                                    
            while n > N or current_state[j] < 0:
                j += 1
                
                if j >= m:
                    return
                #n += -1 +current_state[j-1] + 1
                n += current_state[j-1]
                current_state = current_state[:j-1] + (1, current_state[j]-1) + current_state[j+1:]
                

    def generator_holes(self,m,N):
        if (m == 0):
            yield ()
            return
        current_state = (0,)*m
        n = 0
        while True:
            if (n == N):
                yield current_state
            j = 0
            current_state = (current_state[j]+1,) + current_state[j+1:]
            n += 1
            while n > N or current_state[j] > 1:
                j += 1
                if j >= m:
                    return
                n -= current_state[j-1] - 1
                current_state = current_state[:j-1] + (0, current_state[j]+1) + current_state[j+1:]


        
class fock_space_kron:
    def __init__(self, f1, f2):

        self.f1 = f1
        self.f2 = f2
        
        #1)
        self.modes = f1.modes + f2.modes

        #2)
        self.dimension  = f1.dimension * f2.dimension

        #3)
        self.emptyH = coo_matrix((self.dimension , self.dimension ), dtype = complex).tocsc()

        #4)
        self.eye = sparse.eye(self.dimension ).tocsc()
        
        st1 = f1.statistics
        st2 = f2.statistics
        
        if ((st1 =='Fermi') and (st2 =='Fermi')) or ((st1 =='Fermi_sea') and (st2 =='Fermi')) or ((st1 =='Fermi') and (st2 =='Fermi_sea')):
            self.statistics = 'Fermi'            
            data1 = np.zeros(f1.dimension)
            data2 = np.zeros(f2.dimension)
            
            for i in range(f1.dimension):
                k = sum(f1.occupations(i))
                data1[i] = (-1)**k
            for i in range(f2.dimension):
                k = sum(f2.occupations(i))
                data2[i] = (-1)**k
        #4.1)
            self.parity_f1 = diags(data1)
        #4.2)
            self.parity_f2 = diags(data2)
        
        #5)
            self.annihilate=[]
        
            for k in range(f1.modes):
                self.annihilate.append(sparse.kron(f1.annihilate[k], f2.eye).tocsc())
            for k in range(f1.modes, self.modes):
                self.annihilate.append(sparse.kron(self.parity_f1, f2.annihilate[k-f1.modes]).tocsc())
        #6)
            self.create=[]
        
            for k in range(f1.modes):
                self.create.append(sparse.kron(f1.create[k], f2.eye).tocsc())
            for k in range(f1.modes, self.modes):
                self.create.append(sparse.kron(self.parity_f1, f2.create[k-f1.modes]).tocsc())
                
                
        else:
            self.statistics = 'Non_fermi'
        #5) 
            self.annihilate=[]
        
            for k in range(f1.modes):
                self.annihilate.append(sparse.kron(f1.annihilate[k], f2.eye).tocsc())
            for k in range(f1.modes, self.modes):
                self.annihilate.append(sparse.kron(f1.eye, f2.annihilate[k-f1.modes]).tocsc())
        #6)
            self.create=[]
        
            for k in range(f1.modes):
                self.create.append(sparse.kron(f1.create[k], f2.eye).tocsc())
            for k in range(f1.modes, self.modes):
                self.create.append(sparse.kron(f1.eye, f2.create[k-f1.modes]).tocsc())
        '''    
        
        self.annihilate=[]
        for k in range(f1.modes):
            self.annihilate.append(sparse.kron(f1.annihilate[k], f2.eye).tocsc())
        for k in range(f1.modes, self.modes):
            self.annihilate.append(sparse.kron(f1.eye, f2.annihilate[k-f1.modes]).tocsc())
        
        self.create=[]
        for k in range(f1.modes):
            self.create.append(sparse.kron(f1.create[k], f2.eye).tocsc())
        for k in range(f1.modes, self.modes):
            self.create.append(sparse.kron(f1.eye, f2.create[k-f1.modes]).tocsc())
        '''         
    def sigmax(self, i):
        if (i<self.f1.modes):
             return(sparse.kron(self.f1.sigmax(i), self.f2.eye).tocsc())
        else:
            return(sparse.kron(self.f1.eye, self.f2.sigmax(i-self.f1.modes)).tocsc())
        

        
    def sigmay(self, i):
       
        if (i<self.f1.modes):
             return(sparse.kron(self.f1.sigmay(i), self.f2.eye).tocsc())
        else:
            return(sparse.kron(self.f1.eye, self.f2.sigmay(i-self.f1.modes)).tocsc())

        
    def sigmaz(self, i):
        
        if (i<self.f1.modes):
             return(sparse.kron(self.f1.sigmaz(i), self.f2.eye).tocsc())
        else:
            return(sparse.kron(self.f1.eye, self.f2.sigmaz(i-self.f1.modes)).tocsc())
    
    
    def occupations(self,i):
        
        if i >= self.dimension :
            print('the number is out of range')
        else:
            state = np.zeros(self.modes, dtype = int)
            state[:self.f1.modes] = self.f1.occupations(i//self.f2.dimension)
            state[self.f1.modes:] = self.f2.occupations(i%self.f2.dimension)
            return(state)
        
        
    def index(self,state):
        s = list(state)
        return(self.f1.index(state[:self.f1.modes]) * self.f2.dimension + self.f2.index(state[self.f1.modes:]))
        
def real_time_solver(psi0, dt, tmax, H, Q = None, final_state = None):
    K = psi0.size
    Nt = int(tmax/dt)+1
    psi = [psi0, psi0]
    
    if callable(H) == False:
        H_const = H
        def H(t):
            return(H_const)
     
    if Q == None:
        for i in range(1,Nt):
            psi_iter_old = psi0
            psi_iter = np.zeros(K)
            psi_compare = psi0
            while (LA.norm(psi_iter-psi_compare)>10**(-6)):
                s = psi[1]+psi_iter_old
                psi_iter = psi[1]+ dt/2*(-1j) * H(i*dt+dt/2)@s
                psi_compare = psi_iter_old
                psi_iter_old = psi_iter
            psi[1] = psi_iter
            
        return (psi[1])
        
    elif type(Q) == list:
        l = len(Q)
        results= np.zeros((l,Nt))
        for j in range(l):
            results[j,0] = abs(np.conj(psi0) @ Q[j] @ psi0)
        
        for i in range(1,Nt):
            psi_iter_old = psi0
            psi_iter = np.zeros(K)
            psi_compare = psi0
            while (LA.norm(psi_iter-psi_compare)>10**(-6)):
                s = psi[1]+psi_iter_old
                psi_iter = psi[1]+ dt/2*(-1j) * H(i*dt+dt/2)@s
                psi_compare = psi_iter_old
                psi_iter_old = psi_iter
            psi[1] = psi_iter
            
            for j in range(l):
                results[j,i] = abs(np.conj(psi[1]) @ Q[j] @ psi[1])
        
    else:
        results = np.zeros(Nt)
        results[0] = abs(np.conj(psi0) @ Q @ psi0)
        
        for i in range(1,Nt):
            psi_iter_old = psi0
            psi_iter = np.zeros(K)
            psi_compare = psi0
            while (LA.norm(psi_iter-psi_compare)>10**(-6)):
                s = psi[1]+psi_iter_old
                psi_iter = psi[1]+ dt/2*(-1j) * H(i*dt+dt/2)@s
                psi_compare = psi_iter_old
                psi_iter_old = psi_iter
            psi[1] = psi_iter
            results[i] = abs(np.conj(psi[1]) @ Q @ psi[1])
    if final_state != None :
        return (results, psi[1])
    else:
        return(results)
    
        



        
        

class fermi_sea_joint:
    def __init__(self, max_total_excitations, num_modes1, num_modes2, fermi_sea_modes1, fermi_sea_modes2):
        
        self.statistics = 'Fermi_sea'
        self.max_total_excitations = max_total_excitations #internal parameter
        
        self.modes1 = num_modes1
        self.modes2 = num_modes2
        self.modes = num_modes1 + num_modes2

        self.fermi_sea_modes1 = fermi_sea_modes1
        self.fermi_sea_modes2 = fermi_sea_modes2      
            
        
        self.states_list = list(self.sea_generator())
        
        
        #4)
        self.find_index = {state: idx for idx, state in enumerate(self.states_list)}
        
        #5)
        self.dimension  = len(self.states_list)
        
        #6)
        self.emptyH = coo_matrix((self.dimension , self.dimension ), dtype = complex).tocsc()
        
        #7)
        self.eye = sparse.eye(self.dimension ).tocsc()
        
        
        self.annihilate=[]
        current_state = []
        for k in range(self.modes):
            row = np.zeros(self.dimension )
            col = np.arange(self.dimension , dtype=int)
            data = np.zeros(self.dimension )
            for i in range(self.dimension ):
                if self.states_list[i][k]==0:
                    row[i] = i
                else:
                    current_state = list(self.states_list[i])
                    y = sum(current_state[:k])
                    current_state[k] = 0
                    p = self.index(current_state)
                    
                    if (p=='666'):
                        continue
                    row[i] = p
                    data[i]=(-1)**y*(self.states_list[i][k])**0.5
            A = coo_matrix((data, (row, col)), shape=(self.dimension , self.dimension ), dtype = complex)
            #A.eliminate_zeros()
            self.annihilate.append(A.tocsc())
                              
              
        self.create=[]
        current_state = []
        for k in range(self.modes):
            row = np.zeros(self.dimension )
            col = np.arange(self.dimension , dtype=int)
            data = np.zeros(self.dimension )
            for i in range(self.dimension ):
                if self.states_list[i][k] == 1 :
                    row[i] = i 
                else:
                    current_state = list(self.states_list[i])
                    y = sum(current_state[:k])
                    current_state[k] = 1
                    p = self.index(current_state)
                    if (p=='666'):
                        continue
                    row[i] = p
                    data[i]=(-1)**y*(self.states_list[i][k]+1)**0.5
            A = coo_matrix((data, (row, col)), shape=(self.dimension , self.dimension ), dtype = complex)
            #A.eliminate_zeros()
            self.create.append(A.tocsc())
    
    
    def sigmax(self, i):

        return(self.annihilate[i] + self.create[i])
    
    
    def sigmay(self, i):
 
        return(-1j*(-self.annihilate[i] + self.create[i]))
    
    
    def sigmaz(self, i):

        return(- self.annihilate[i]@self.create[i] + self.create[i]@self.annihilate[i])
    
    
    def occupations(self,i):
     
        if i >= self.dimension :
            print('the number is out of range')
            
        else:
            return(np.array(self.states_list[i]))
        
    
    
    def index(self, state):
       
        if len(state) != self.modes:
            print('incorrect size of an array')
        else:
            s = tuple(state)
            try:
                return(self.find_index[s])
            except:
                return('666')
  
        
            
    def occupations(self,i):
     
        if i >= self.dimension :
            print('the number is out of range')
            
        else:
            return(np.array(self.states_list[i]))
    
    
            
       
    
    def generator_particles(self,m,N):
        if (m == 0):
            yield ()
            return
        current_state = (1,)*m
        n = 0
        while True:
            if (n == N):
                yield current_state
            j = 0
            current_state = (current_state[j]-1,) + current_state[j+1:]
            n += 1
                                    
            while n > N or current_state[j] < 0:
                j += 1
                
                if j >= m:
                    return
                #n += -1 +current_state[j-1] + 1
                n += current_state[j-1]
                current_state = current_state[:j-1] + (1, current_state[j]-1) + current_state[j+1:]
                

    def generator_holes(self,m,N):
        if (m == 0):
            yield ()
            return
        current_state = (0,)*m
        n = 0
        while True:
            if (n == N):
                yield current_state
            j = 0
            current_state = (current_state[j]+1,) + current_state[j+1:]
            n += 1
            while n > N or current_state[j] > 1:
                j += 1
                if j >= m:
                    return
                n -= current_state[j-1] - 1
                current_state = current_state[:j-1] + (0, current_state[j]+1) + current_state[j+1:]


    def sea_generator(self):
        m_l = self.modes1
        mp_l = self.fermi_sea_modes1
        mh_l = m_l - mp_l
        m_r = self.modes2
        mp_r = self.fermi_sea_modes2
        mh_r = m_r - mp_r
        
        N = self.max_total_excitations
        
        
        
        if N>(m_l+m_r):
            N = m_l+m_r
        for Nl in range(N+1):
            Fl1 = min(Nl,mp_l)
            Fl2 = min(Nl,mh_l)
            for Nr in range(N-min(Nl,m_l)+1):
                Fr1 = min(Nr,mp_r)
                Fr2 = min(Nr,mh_r)
                
                for i in range (Nl-Fl1,Fl2+1):
                    for j in range (Nr-Fr1,Fr2+1):
                        for hl1 in self.generator_particles(mp_l, Nl-i):    
                            for hl2 in self.generator_holes(mh_l,i):
                                for hr1 in self.generator_particles(mp_r,Nr-j):
                                    for hr2 in self.generator_holes(mh_r,j):
                                        yield (hl1 + hl2 + hr1 + hr2)

                
                    
                    

                            
def real_time_solver1(psi0, psi_vac, dt, tmax, H, Q = None, final_state = None):
    K = psi0.size
    Nt = int(tmax/dt)+1
    psi = [psi0, psi0]
    
    if callable(H) == False:
        H_const = H
        def H(t):
            return(H_const)
     
    if Q == None:
        psi_ar = np.zeros((K,Nt), dtype = complex)
        psi_ar[:,0] = psi0
        for i in range(1,Nt):
            psi_iter_old = psi0
            psi_iter = np.zeros(K)
            psi_compare = psi0
            while (LA.norm(psi_iter-psi_compare)>10**(-6)):
                s = psi[1]+psi_iter_old
                psi_iter = psi[1]+ dt/2*(-1j) * H(i*dt+dt/2)@s
                psi_compare = psi_iter_old
                psi_iter_old = psi_iter
            psi[1] = psi_iter
            psi_ar[:,i] = psi[1]
            
        return (psi_ar)
        
    elif type(Q) == list:
        l = len(Q)
        results= np.zeros((l,Nt))
        for j in range(l):
            #results[j,0] = np.conj(psi_vac) @ Q[j] @ psi0
            results[j,0] = (Q[j] @ psi0)[0]
        
        for i in range(1,Nt):
            psi_iter_old = psi0
            psi_iter = np.zeros(K)
            psi_compare = psi0
            while (LA.norm(psi_iter-psi_compare)>10**(-6)):
                s = psi[1]+psi_iter_old
                psi_iter = psi[1]+ dt/2*(-1j) * H(i*dt+dt/2)@s
                psi_compare = psi_iter_old
                psi_iter_old = psi_iter
            psi[1] = psi_iter
            
            for j in range(l):
                #results[j,i] = np.conj(psi_vac) @ Q[j] @ psi[1]
                results[j,i] = (Q[j] @ psi[1])[0]
        
    else:
        results = np.zeros(Nt)
        results[0] = abs(np.conj(psi_vac) @ Q @ psi0)
        
        for i in range(1,Nt):
            psi_iter_old = psi0
            psi_iter = np.zeros(K)
            psi_compare = psi0
            while (LA.norm(psi_iter-psi_compare)>10**(-6)):
                s = psi[1]+psi_iter_old
                psi_iter = psi[1]+ dt/2*(-1j) * H(i*dt+dt/2)@s
                psi_compare = psi_iter_old
                psi_iter_old = psi_iter
            psi[1] = psi_iter
            results[i] = abs(np.conj(psi[1]) @ Q @ psi[1])
    if final_state != None :
        return (results, psi[1])
    else:
        return(results)    