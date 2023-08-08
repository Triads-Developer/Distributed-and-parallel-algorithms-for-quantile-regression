import qpadmslack
import time
import random
import concurrent.futures
import numpy as np
from scipy.sparse import linalg as sla

def main():
    N = 100000
    p = 100
    tau = 0.7
    rho = 0.5
    a = 3.7
    pho = 5
    beta_true = [0] * p
    beta_true[6] = beta_true[12] = beta_true[15] = beta_true[20] = 1
    #hardcoding beta_true[1]
    beta_true[1] = 0.367
    #beta_true[1] = 0.7*qnorm(tau)
    ##we consider five different values for the partition number K, i.e., 1, 10, 20, 50 and 100
    K = [1, 10, 20, 50, 100]
    ##the corresponding suitable value for lambda under different K
    lamb = [20, 30, 60, 100, 200]
    ##AE and Time respectively record the estimation accuracy and computational time (pseudo-parallel) of the algorithm under different K
    testxk = AE = Time = [0] * len(K)

    start = time.perf_counter()
    fullthreadedPara = qpadmslack.fullthreadedparaQPADMslackcpp([1.2,3.2,0.2345],100,tau,"scad",a,200,pho,1000,0.001,False)
    end = time.perf_counter()
    print(f'The threaded Python call with a large input finished in {round(end-start, 2)} second(s) with a K value of 100 ')
    #This is loading the file in instead of generating the data
    #X = np.array(np.loadtxt("../data/X.backup", delimiter=",", dtype="float32"))
   
    #for k in range(len(K)):
    #    start = time.perf_counter()
    #    fullthreadedPara = qpadmslack.fullthreadedparaQPADMslackcpp([1.2,3.2,0.2345],K[k],tau,"scad",a,lamb[k],pho,1000,0.001,False)
    #    end = time.perf_counter()
    #    print(f'The threaded Python call with a large input finished in {round(end-start, 2)} second(s) with a K value of {K[k]}')
  
    #para= qpadmslack.paraQPADMslackcpp(X,K[2],tau,"scad",a,lamb[2],pho,1000,0.001,False)
    #fullthreadedPara = qpadmslack.fullthreadedparaQPADMslackcpp(X,K[2],tau,"scad",a,lamb[2],pho,1000,0.001,False)
    #print("fullthreaded version " , fullthreadedPara)
if __name__ == '__main__':
    main()
