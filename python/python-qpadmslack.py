import qpadmslack
import sys, getopt
import time
import random
import concurrent.futures
import numpy as np
from scipy.sparse import linalg as sla

def printHelp():
    print ('python-qpadmslack.py \n-x <path_to_input_file_x> \n-y <path_to_input_file_x> \n-t (tau) \n-r (rho) \n-a \n-p (pho) \n-k \n-l (lamb) \n-m (maxstep) \n-e (eps) \n-i (intercept) \n-s (saveBetaHistory) \n-p (penalty)')

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hx:y:t:r:a:p:k:l:m:e:i:s:p",[])
        if len(opts) == 0:
            printHelp()
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            printHelp()
            sys.exit()
        elif opt in ("-x", "--xfile"):
            pathToX = arg
        elif opt in ("-y", "--yfile"):
            pathToY = arg
        elif opt in ("-tau"):
            tau = float(arg)
        elif opt in ("-rho"):
            rho = float(arg)
        elif opt in ("-a"):
            a = float(arg)
        elif opt in ("-pho"):
            pho = float(arg)
        elif opt in ("-k"):
            k = int(arg)
        elif opt in ("-lamb"):
            lamb = int(arg)
        elif opt in ("-maxstep"):
            maxstep = int(arg)
        elif opt in ("-eps"):
            eps = float(arg)
        elif opt in ("-intercept"):
            intercept = arg
        elif opt in ("-saveBetaHistory"):
            saveBetaHistory = arg == "True"
        elif opt in ("-penalty"):
            penalty = arg == "True"

   #this is just a default value for X and Y's path - it's more just
   #convencience for me and testing
    try: pathToX
    except NameError: pathToX = "../data/X.small"

    try: pathToY
    except NameError: pathToY = "../data/Y.small"

    try: tau
    except NameError: tau = 0.7

    try: penalty
    except NameError: penalty = "scad"

    try: rho
    except NameError: rho = 0.5

    try: a
    except NameError: a = 3.7

    try: pho
    except NameError: pho = 5

    ##we consider five different values for the partition number K, i.e., 1, 10, 20, 50 and 100
    #K = [1, 10, 20, 50, 100]
    ##the corresponding suitable value for lambda under different K
    #lamb = [20, 30, 60, 100, 200]
    try: k
    except NameError: k = 100

    try: lamb
    except NameError: lamb = 200

    try: maxstep
    except NameError: maxstep = 1000

    try: eps
    except NameError: eps= 0.001
    
    try: intercept
    except NameError: intercept= False

    try: saveBetaHistory
    except NameError: saveBetaHistory= False
    
    fullthreadedPara = qpadmslack.fullthreadedparaQPADMslackcpp(pathToX, pathToY, k,tau,penalty,a,lamb,pho,maxstep,eps,intercept,saveBetaHistory)

    print(f'{fullthreadedPara[0][0,0]} second(s)\n')
    print(f'{fullthreadedPara[0][0,1]} second(s)\n')
    print(f'{fullthreadedPara[0][0,2]} second(s) to load the data\n')
    print(f'{fullthreadedPara[0][0,3]} second(s) to perform the matrix inversions \n')
    print(f'{fullthreadedPara[0][0,4]} second(s) to perform the iterative portion of the algorithm\n')
    print(f'{fullthreadedPara[1]} beta\n ')
    print(f'{fullthreadedPara[2]} iterations\n ')
    if saveBetaHistory == True:
        print(f'{fullthreadedPara[3]} beta values for each iteration\n ')

if __name__ == '__main__':
    main(sys.argv[1:])
