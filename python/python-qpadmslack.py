import qpadmslack
import sys, getopt
import time
import random
import concurrent.futures
import numpy as np
from scipy.sparse import linalg as sla

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hx:y:",["xfile=, yfile="])
        if len(opts) == 0:
            print ('python-qpadmslack.py -x <inputXfile> -y <inputY>')
            sys.exit(2)
    except getopt.GetoptError:
        print ('python-qpadmslack.py -x <inputXfile> -y <inputY>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print ('python-qpadmslack.py -x <inputXfile> -y <inputY>')
            sys.exit()
        elif opt in ("-x", "--xfile"):
            pathToX = arg
        elif opt in ("-y", "--yfile"):
            pathToY = arg

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

    fullthreadedPara = qpadmslack.fullthreadedparaQPADMslackcpp(pathToX, pathToY, [1.2,3.2,0.2345],100,tau,"scad",a,200,pho,1000,0.001,False)
    print(f'The threaded Python call with an input finished in {fullthreadedPara[0]} second(s) with a K value of 100 ')

if __name__ == '__main__':
    main(sys.argv[1:])
