# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 16:10:58 2021

@author: kst

Toy example for garbled circuit for evaluation of max(x-y)
"""


import numpy as np
import random
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
import time
from threading import Thread
import queue as que

from GarblingFreeXOR_c import garble
from Elgamal import ElGamal
from obliviousTransfer import OT
from FFArithmetic import GF

WIRE_BITLENGTH = 128
class cloud(Thread):
    r'''
    Implement the functionality of each of the two clouds. 
    '''
    def __init__(self, L, M, b, c, x, triplets, elgamal, inQue, outQue, bitlen, num):
        Thread.__init__(self)
        self.elg = elgamal                  #used for oblivious transfer
        self.L = L
        self.M = M
        self.b = b
        self.c = c
        self.x = x
        self.triplets = triplets
        
        self.inQue = inQue
        self.outQue = outQue
        self.bitlen = bitlen
        self.num = num   # Number of the cloud (either 0 or 1)
        
        self.GARB = garble(WIRE_BITLENGTH, self.bitlen)
        self.OT = OT(bitlen, elgamal, p1)
        
    def SSmultParal(self,A,b,triplets):
        r'''
        Multiplikation of shares using Beavers triplets. 
        '''
        n,m = np.shape(A)
        out = []
        A1 = A.flatten()
    
        for i in range(n*m):
            out.append( [ (A1[i] - triplets[i][0]) %p1 , (b[i % m] - triplets[i][1]) %p1 ] )
        
        self.outQue.put(out)
        
        rec = self.inQue.get(block = True)
        
        de = [((i[0] + j[0]) % p1, (i[1]+j[1]) % p1) for i,j in zip(out,rec)]
    
        A2 = []
        deProd = []
        for i,tub in enumerate(de):
            A2.append( (  int(tub[0]) * int(triplets[i][1]) + int(tub[1])* int(triplets[i][0])  + int(triplets[i][2])) % p1)
            deProd.append(tub[0]*tub[1] % p1)
            
        A2 = np.array(A2).reshape((n,m))
        
        res = np.sum(A2,axis=1) % p1
        
        res2 =( res + np.sum(np.array(deProd).reshape((n,m)), axis = 1)) % p1
        
        return res, res2
        
    def convertInp(self, inp):
        r'''
        Convert int to bin and save each list of binary input in a dictionary.
        '''
        inputs = {}
        for i in range(len(inp)):                                      
            temp = []
            for j in range(self.bitlen):
                temp.append(int(bin(inp[i])[2:].rjust(self.bitlen,'0')[j]) )
            inputs[i] = temp[::-1]
        return inputs

    def run(self):
        
        R =random.getrandbits(self.bitlen) 
#        if self.num == 0:
#            print('R', R)
        n1,m1 = np.shape(self.L)
        n2,m2 = np.shape(self.M)
        n = [n1,n2]
        # Calculate Lx + b and Mx + c using additive secret sharing
        
        Lx   = self.SSmultParal(self.L, self.x, self.triplets[:n1*m1])[self.num]
        Lx_b = (Lx.reshape(n1,1) + self.b) % p1
        
        Mx   = self.SSmultParal(self.M, self.x, self.triplets[n1*m1:])[self.num]
        Mx_c = (Mx.reshape(n2,1) + self.c) % p1
        
        
        #Calculate inputs to garbling
        inp1 = Lx_b.reshape((n1,))
        inp2 = Mx_c.reshape((n2,))
        
        if self.num == 0:
            inp2 = np.append(inp2, R)
        elif self.num == 1:
            inp1 = np.append(inp1, R)

        # Convert the inputs from int to bin and arrange the inputs in a dictionary
        inputs1 = self.convertInp(inp1)
        inputs2 = self.convertInp(inp2)
        #Use the list so each cloud can get the correct input for garbling (since they do one garbling each in parallel)
        inputslist = [inputs1, inputs2]

        pk = self.OT.choose(inputslist[self.num])
        self.outQue.put(pk)
        # Start the garbling of one of the inputs because the other cloud does the other input.
#        if self.num == 0:
        F, Y, d, LEFTinp = self.GARB.garble(inputslist[(self.num+1)%2])
        self.F = F
        pk = self.inQue.get(block = True)
        c = self.OT.transfer(pk, LEFTinp)
        self.outQue.put(c)
        self.outQue.put([F,Y,d])
        
#        if self.num == 1:
#            pk = self.inQue.get(block = True)
            #Evaluate the received garbled circuit
        c = self.inQue.get(block = True)
        X = self.OT.retrieve(c, inputslist[self.num])
        F,Y,d = self.inQue.get(block = True)
        res = self.GARB.evaluate(F,Y,X,d, n[self.num])
        res = binToint(self.bitlen, res)
#        if self.num == 0:
#            print(res)
        self.res =  TwoCompl(self.bitlen, res + R)
#        if self.num == 0:
#            print('res',  self.res)
            
#        if self.num == 1:
#            time.sleep(1)
##            print(inputslist)
#            print('R', R)
#            print('res', res)
#            print(self.res)
        
def TwoCompl(bitlen, intnum):
    if intnum > 2**(bitlen-1):
        return intnum - 2**bitlen
    elif intnum < -2**(bitlen-1):
        return intnum + 2**bitlen
    else:
        return intnum

def binToint(bitlen, binnum):
    '''
    Convert bin to int and also decide if number is negative
    '''
    r = ''
    for i in binnum:
        r+=str(i)
    intnum = TwoCompl(bitlen, int(r,2))
    return intnum

                


def share(A):
    '''
    Generate two additive shares of each element in A
    '''
    n,m = np.shape(A)
    
    A1 = np.random.randint(0,p1,size=(n,m), dtype = np.int64)
    A2 = (A - A1) % p1
    
    return A1,A2

def gen_triplets(n):
    '''
    Generate n triplets and arrange them in a list where the first element are one of the shares of all triplets and 
    second element the other share of all triplets.
    '''
    
    a = np.random.randint(0,p1,(n,1),dtype=np.int64)
    a1,a2 = share(a)
    b = np.random.randint(0,p1, (n,1),dtype=np.int64)
    b1,b2 = share(b)
    c = (a*b) % p1
    c1,c2 = share(c)
    return [(i[0],j[0],l[0]) for i,j,l in zip(a1,b1,c1)]  ,  [(i[0],j[0],l[0]) for i,j,l in zip(a2,b2,c2)] 
    

#################### MAIN STUFF ##################################
#sss = 4
#np.random.seed(sss)
#random.seed(sss)
#bitlen = 5                                             # SHOULD NOT BE LARGER THAN 16. input bitlen. inputs go from -2**(bitlen-1) : 2**(bitlen -1)-1 
#bitlenRes = 2*bitlen                                    #bitlen after multiplication and addition.
#WIRE_BITLENGTH = 128                                    #privacy parameter for the garbling

## Secret sharing parameters
#p1 = 1967753     #21 bit prime
#p1 = 2**(bitlenRes) #690444331  # 30 bit prime          #For secret sharing we do modular the number

## Elgamal parameters
p = 101518638042068819827067989441957810604363          #A huge prime for Elgamal
G = GF(p)                                               #Group with DDH asumption, has to be huge because the wire labels are 128 bit.
q = 50759319021034409913533994720978905302181           #Another huge prime 
g = G(2)                                                #Generating element for the group.
elg   = ElGamal(G,g,q)                                  #Instance of Elgamal class        

## Generation of problem matrices
GAT = []
#n3 = np.array([2,4,6,8,10,12,14])#range(3,11)
#X,Y = np.meshgrid(n3,n3)
#for iii in n3:
#    for jjj in n3:
l = 16
bitlen = l/2                                           # SHOULD NOT BE LARGER THAN 16. input bitlen. inputs go from -2**(bitlen-1) : 2**(bitlen -1)-1 
bitlenRes = l   
p1 = 2**(l) 

nn = 16 #needs to be power of 2
n1,m = (nn,2)                                             # Dimension of the matrices     
n2,m = (nn,2)
L = np.random.randint(-2**(bitlen-1)+1, 2**(bitlen-1), size=(n1,m))
M = np.random.randint(-2**(bitlen-1)+1, 2**(bitlen-1), size=(n2,m))
b = np.random.randint(-2**(bitlen-1)+1, 2**(bitlen-1), size=(n1,1))
c = np.random.randint(-2**(bitlen-1)+1, 2**(bitlen-1), size=(n2,1))
x = np.random.randint(1, 2**(bitlen-1), size=(m,1))

#generate 2 shares of each parameter
L1,L2 = share(L)
M1,M2 = share(M)
b1,b2 = share(b)
c1,c2 = share(c)
x1,x2 = share(x)
#print(np.dot(L,x)+b)
#print( np.dot(M,x)+c )
#generate a triple for each multiplication
t1, t2 = gen_triplets(m*(n1+n2))

#Some queues for communication bewteen the clouds
Que1 = que.Queue()
Que2 = que.Queue()

#Instanciate two clouds, with number 0 and 1 and each their own shares of the parameters
cloud1 = cloud(L1,M1,b1,c1,x1,t1, elg, Que1, Que2, bitlenRes, 0)
cloud2 = cloud(L2,M2,b2,c2,x2,t2, elg, Que2, Que1, bitlenRes, 1)

#Start the computation and time it
time1 = time.time()
cloud1.start()
cloud2.start()

cloud1.join()
cloud2.join()
time2 = time.time()

F = cloud1.F

#Get the results from the threads (clouds)
res1 = cloud1.res
res2 = cloud2.res

#Get the garbled circuit just to see the size of it.
#F = Que1.get(block = True)

gateNum = 0
for i in F[:-(nn+l)+1]:
    gateNum += 1#len(i) % 3

for i in F[-(nn+l)+1:-l]:
    for j in i:
        gateNum += len(j) % 3
for i in F[-l:]:
    gateNum += 1#len(i) % 3


###Convert the result from bin to int and compute the subtraction
res = TwoCompl(l, res1 -res2)

#Print the results

print('Expected result: ', int(max(np.dot(L,x)+b) - max(np.dot(M,x)+c)  )   )
print('Garbled result:  ', res)
print('Execution time:  ', time2-time1, 'seconds')
print('Size of GC:      ', str(sys.getsizeof(F)), 'bytes')
print('# of gates in GC:', gateNum, 'gates')
print('Calc # of gates: ', 4*nn*l -2*l, 'gates')
print('Exe time / gate: ', ((time2-time1) / gateNum)*1000, 'ms')
##fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot_surface(X, Y, GAT,
#                       linewidth=0, antialiased=False)
#plt.plot(n3,GAT)