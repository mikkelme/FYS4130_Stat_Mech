import numpy as np 
import matplotlib.pyplot as plt
from plot_set import *
from matplotlib.ticker import MaxNLocator
import subprocess

class analyze:
    def  __init__(self, filename, fetch_file = False):
        self.filename = filename
        if fetch_file: fetch_file_from_server(filename)
        self.read_data()

    def read_data(self):
        infile = open(self.filename, 'r')

        info = []
        for i in range(7):
            info.append(int(infile.readline().split()[-1].strip('\n')))
        self.dim, self.L, self.q, self.cluster_per_cycle, self.MC_equil, self.MC_measure, self.bins = info
         
        T_tmp = np.zeros(self.bins)
        m_tmp = np.zeros(self.bins)
        m1_tmp = np.zeros(self.bins)
        m2_tmp = np.zeros(self.bins)
        m4_tmp = np.zeros(self.bins)
        C_r_tmp = np.zeros((self.bins, self.L**self.dim))
        gamma_tmp = np.zeros(self.bins)

        T = []; m = []; m1 = []; m2 = []; m4 = []; C_r = []; gamma = []
        std_T = []; std_m = []; std_m1 = []; std_m2 = []; std_m4 = []; std_C_r = []; std_gamma = []
        
        next(infile) # skip header line
        line_count = 0
        for line in infile: 
            i = line_count % self.bins
            words = line.split()
            T_tmp[i] = float(words[0]) 
            m_tmp[i] = float(words[1])
            m1_tmp[i] = float(words[2]) 
            m2_tmp[i] = float(words[3]) 
            m4_tmp[i] = float(words[4]) 
            C_r_tmp[i] = words[5:]
            gamma_tmp[i] = m4_tmp[i]/m2_tmp[i]**2

            if i == self.bins - 1:
                # Mean value 
                T.append(np.mean(T_tmp))
                m.append(np.mean(m_tmp))
                m1.append(np.mean(m1_tmp))
                m2.append(np.mean(m2_tmp))
                m4.append(np.mean(m4_tmp))
                C_r.append(np.mean(C_r_tmp, axis = 0))
                gamma.append(np.mean(gamma_tmp))

                # Standard deviation of the mean 
                std_T.append(np.std(T_tmp, ddof = 1))
                std_m.append(np.std(m_tmp, ddof = 1))
                std_m1.append(np.std(m1_tmp, ddof = 1))
                std_m2.append(np.std(m2_tmp, ddof = 1))
                std_m4.append(np.std(m4_tmp, ddof = 1))
                std_C_r.append(np.std(C_r_tmp, axis = 0, ddof = 1))
                std_gamma.append(np.std(gamma_tmp, ddof = 1))

            line_count += 1

        self.T = np.array(T)
        self.m = np.array(m)
        self.m1 = np.array(m1)
        self.m2 = np.array(m2)
        self.m4 = np.array(m4)
        self.C_r = np.array(C_r, dtype = float)
        self.gamma = self.m4 / self.m2**2 

        self.std_T = np.array(std_T)
        self.std_m = np.array(std_m)
        self.std_m1 = np.array(std_m1)
        self.std_m2 = np.array(std_m2)
        self.std_m4 = np.array(std_m4)
        self.std_C_r = np.array(std_C_r, dtype = float)
        self.std_gamma = np.array(std_gamma)



def fetch_file_from_server(filenames):   
    ssh = 'mikkelme@ml6.hpc.uio.no'
    remote_path = '/itf-fi-ml/home/mikkelme/FYS4130/'
    destination = './'
    if isinstance(filenames, list):
        filenames = "{" + str(filenames).strip('[]').replace(' ', '').replace('\'', '') + "}"
     
    fetch = input(f'Download: \"{filenames}\" from ssh (y / RN key): ')
    if fetch == 'y': 
        remote = ssh + ":" + remote_path + filenames
        subprocess.run(["scp", remote, destination])
    else:
        print('Aborted download')


def C(r, T, L):
    a = (np.exp(1/T) - 1)
    b = (np.exp(1/T) + 2)
    nominator = a**L + a**r * b**(L-r) + a**(L-r) * b**r
    denominator = 2*a**L + b**L
    return nominator/denominator

if __name__ == '__main__':
    fetch_file_from_server(['2D_chain_2f_L8.txt','2D_chain_2f_L16.txt', '2D_chain_2f_L32.txt'])

  
