# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 16:10:58 2021

@author: kst

Toy example for garbled circuit for evaluation of max(x-y)
"""


import numpy as np
import random
import matplotlib.pyplot as plt
import sys
import time
from threading import Thread
import queue as que

#from Garbling_c import garble
from GarblingFreeXOR_c import garble
from Elgamal import ElGamal
from obliviousTransfer import OT
from FFArithmetic import GF


# controller parameters for 8 segments (either choose 8 or 16)
#N=8;
#L= np.array([[-0.0738 ,0.3096, 0.0128,-0.0065,0.3096,0.0088 ,0.0688,-0.3120],[ -0.5160,-0.3165,-0.2996, -0.3997,  0.6835,  0.5231, -0.4057,-0.6389]]).transpose()
#b= np.array( [-0.3713, -4.6018 ,-0.4953, -0.5374, -0.6019, -1.6653, -1.1442, 0.3944]).reshape(N,1)
#M= np.array([[0.3096,0.1179,-0.0738,-0.3120,-0.3120, -0.0804 ,  0.0088,0.0807],[0.6835, -0.0303,-0.5160, 0.3611,-0.6389, 0.4768,0.5231,0.0127]]).transpose()
#c=np.array([0.3981, -0.8778,-1.3712,-4.6055, -0.6056, -1.3885,   -0.6653,   -0.4034]).reshape(N,1)

# controller parameters for 16 segments
N=16;
L=np.array([[ -0.3430,-0.2671,-0.0117, 0.2079 ,0.2413 ,0.3311 ,-0.3275,  0.3311,0.0055 ,-0.0628 ,  -0.3291,  -0.1891 ,  0.0296 ,0.1010 ,0.0316 ,-0.1378],\
            [ -0.7938,-0.4726,-0.9342,0.1099,0.5993, 0.5030,-0.4526,-0.4970, -0.6671,0.3074,-0.8219, -0.8143,0.3518,-0.6149, 0.3579,-0.7303]]).transpose()
b=np.array([0.7014, 0.4947,-0.8832,-0.8494,-0.1865,0.3658,-0.4385,-3.6341, 0.3215,-1.0671,1.3659, 0.9310,-0.4805,  -0.4469,-0.4910, 0.9235]).reshape(N,1)
M=np.array([[0.0316, -0.1891, -0.3291 ,-0.0628 , -0.2145,0.3311, -0.3411 ,-0.0222 ,0.2413,0.0055, -0.3291 ,0.1398 , 0.3438  , 0.1285,-0.2722 ,-0.1178],\
            [0.3580,  -0.8143, 0.1781,0.3079,-0.2734, 0.5030,-0.6308,-0.4318,0.5993,-0.6671,-0.8219, 0.4114,0.5083,-0.1253,  -0.5634,  0.2835]]).transpose()
c =np.array([0.5089,-0.0690,-3.6340,-0.0692,-0.1093, 1.3658,-0.8367,-0.6831,0.8135,-0.6785, 0.3659, 0.9235, 1.0489, 0.1732,-0.9257,-0.6730]).reshape(N,1)

class cloud(Thread):
    r'''
    Implement the functionality of each of the two clouds. 
    '''
    def __init__(self, L, M, b, c, x, triplets, elgamal, inQue, outQue, senQue, bitlen, num, ite):
        Thread.__init__(self)
        self.elg = elgamal                  #used for oblivious transfer
        self.L = L
        self.M = M
        self.b = b
        self.c = c
        self.x = x
        self.triplets = triplets
        self.ite = ite
        
        self.inQue = inQue
        self.outQue = outQue
        self.senQue = senQue
        self.bitlen = bitlen
        self.num = num   # Number of the cloud (either 0 or 1)
        
        self.GARB = garble(WIRE_BITLENGTH, self.bitlen)
        self.OT = OT(bitlen, elgamal, p1)
        
    def SSmultParal(self,A,b,triplets):
        r'''
        Multiplikation of shares using Beavers triplets. 
        '''
        n,m = np.shape(A)
        B = np.repeat(b,n).reshape(m,n).transpose()

        out1 = A - triplets[0]
        out2 = B - triplets[1]

        self.outQue.put([out1,out2])
        
        rec = self.inQue.get(block = True)
        d = out1 + rec[0]
        e = out2 + rec[1]

        res = (np.multiply(d, triplets[1]) + np.multiply(e,triplets[0]) + triplets[2]) % p1

        res1 = np.dot(res.astype(int), np.ones((m,1),dtype=int)) % p1
        res2 = (res + d*e) % p1
        res2 = np.dot(res2, np.ones((m,1),dtype = int)) % p1

        
        return res1, res2
        
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
        
        n1,m1 = np.shape(self.L)
        n2,m2 = np.shape(self.M)
        n = [n1,n2]
        
        time1 = []
        time2 = []
        # Calculate Lx + b and Mx + c using additive secret sharing
        for i in range(self.ite):
            #R is used to 
            R = random.getrandbits(self.bitlen) 
            if i != 0:
                self.x = self.inQue.get(block = True)
            time1.append(time.time())
            Lx   = self.SSmultParal(self.L, self.x, self.triplets[0])[self.num]
            Lx_b = (Lx.reshape(n1,1) + self.b) % p1
            
            Mx   = self.SSmultParal(self.M, self.x, self.triplets[1])[self.num]
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
            F, Y, d, LEFTinp = self.GARB.garble(inputslist[(self.num+1)%2])
            pk = self.inQue.get(block = True)
            c = self.OT.transfer(pk, LEFTinp)
            self.outQue.put(c)
            self.outQue.put([F,Y,d])
            
            #Evaluate the received garbled circuit
            c = self.inQue.get(block = True)
            X = self.OT.retrieve(c, inputslist[self.num])
            F,Y,d = self.inQue.get(block = True)
            res = self.GARB.evaluate(F,Y,X,d, n[self.num])
            res = binToint(self.bitlen, res)
            self.res =  TwoCompl(self.bitlen, res + R)
        
            self.senQue.put( self.res )
            time2.append(time.time())
        self.outQue.put( np.average(np.array(time2) - np.array(time1)) )
        self.F = F

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

def gen_triplets(n,m):
    '''
    Generate n triplets and arrange them in a list where the first element are one of the shares of all triplets and 
    second element the other share of all triplets.
    '''
    
    a = np.random.randint(0,p1,(n,m),dtype=np.int64)
    a1,a2 = share(a)
    b = np.random.randint(0,p1, (n,m),dtype=np.int64)
    b1,b2 = share(b)
    c = np.multiply(a,b) % p1
    c1,c2 = share(c)
    return [a1,b1,c1], [a2,b2,c2] # [(i[0],j[0],l[0]) for i,j,l in zip(a1,b1,c1)]  ,  [(i[0],j[0],l[0]) for i,j,l in zip(a2,b2,c2)] 
    

def f(A,B,x, u):
    return np.dot(A,x) + B*u
    
#################### MAIN STUFF ##################################

#np.random.seed(2)
l=32 # 16 or 32 bit
bitlen = 16                                             # SHOULD NOT BE LARGER THAN 16. input bitlen. inputs go from -2**(bitlen-1) : 2**(bitlen -1)-1 
bitlenRes = 2*bitlen                                    #bitlen after multiplication and addition.
WIRE_BITLENGTH = 128                                    #privacy parameter for the garbling
ite = 21                                                 #number of "iterations" for the controller        

## Secret sharing parameters
#p1 = 1967753     #21 bit prime
p1 = 2**(l) #690444331  # 30 bit prime          #For secret sharing we do modular the number

## Elgamal parameters
p = 101518638042068819827067989441957810604363          #A huge prime for Elgamal
G = GF(p)                                               #Group with DDH asumption, has to be huge because the wire labels are 128 bit.
q = 50759319021034409913533994720978905302181           #Another huge prime 
g = G(2)                                                #Generating element for the group.
elg   = ElGamal(G,g,q)                                  #Instance of Elgamal class        

## Generation of problem matrices
#n1,m = (8,2)                                             # Dimension of the matrices     
#n2,m = (8,2)
#L = np.random.randint(-2**(bitlen-1), 2**(bitlen-1), size=(n1,m))
#M = np.random.randint(-2**(bitlen-1), 2**(bitlen-1), size=(n2,m))
#b = np.random.randint(-2**(bitlen-1), 2**(bitlen-1), size=(n1,1))
#c = np.random.randint(-2**(bitlen-1), 2**(bitlen-1), size=(n2,1))
#x = np.random.randint(0, 2**(bitlen-1), size=(m,1))


bit=(2**(l-1)-1)
biggestSegmentG=13.9251 # from linear programs (not here)
biggestSegmentH=14.2114 
biggestNum=max([biggestSegmentG,biggestSegmentH]) # select max number

s0s1=bit/biggestNum;                # big s for accuracy
s0= int(np.sqrt(s0s1));               # choice: s0=s1;
s1=s0;
s2=s0*s1;          
## Example from paper:
n,m = np.shape(L)
A = np.array([[1,1],[0,1]])
B = np.array([0.5,1]).reshape(2,1)
L = s1 * L
M = s1 * M
b = s2 * b
c = s2 * c
x = np.array([-22,1]).reshape(m,1)

#Rounding
L = L.astype(int)
M = M.astype(int)
b = b.astype(int)
c = c.astype(int)

#generate 2 shares of each parameter
L1,L2 = share(L)
M1,M2 = share(M)
b1,b2 = share(b)
c1,c2 = share(c)
x1,x2 = share(s0*x)

#generate a triple for each multiplication
t1, t2 = gen_triplets(n,m)
t11, t22 = gen_triplets(n,m)
t = [[t1,t11], [t2,t22]]
#Some queues for communication bewteen the clouds
Que1 = que.Queue()
Que2 = que.Queue()
QueA = que.Queue()
QueB = que.Queue()
#Instanciate two clouds, with number 0 and 1 and each their own shares of the parameters
cloud1 = cloud(L1,M1,b1,c1,x1,t[0], elg, Que1, Que2, QueA, l, 0, ite)
cloud2 = cloud(L2,M2,b2,c2,x2,t[1], elg, Que2, Que1, QueB, l, 1, ite)
#
##Start the computation and time it
time1 = []
time2 = []
cloud1.start()
cloud2.start()

u = []
xk = np.zeros((m,ite+1))
xk[:,0] = x.reshape(m,)
time1.append(time.time()) 
for i in range(ite):
#    print(xk[i])
    res1 = QueA.get(block = True)
    res2 = QueB.get(block = True)
    u.append(  TwoCompl(l, res1 - res2) / s2 )
    
    x_k = ( f(A,B,xk[:,i].reshape(m,1),u[i]) )
    xk[:,i+1] = x_k.reshape(m,)
    x_k*=s0
    x1,x2 = share(x_k.astype(int))
    
    Que1.put(x1)
    Que2.put(x2)
   

cloud1.join()
cloud2.join()
time2.append(time.time())

TIME = np.average(np.array(time2)-np.array(time1))

plt.plot(xk[0,:], xk[1,:])
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.show()

#Get the results from the threads (clouds)
TIME2 = Que1.get(block = True)
TIME3 = Que2.get(block = True)

print('Execution: ', TIME )
print('gen of u(k): ', 1/2 * (TIME2 + TIME3))

F = cloud1.F

gateNum = 0
for i in F[:-(n+l)+1]:
    gateNum += 1#len(i) % 3

for i in F[-(n+l)+1:-l]:
    for j in i:
        gateNum += len(j) % 3
for i in F[-l:]:
    gateNum += 1#len(i) % 3

print('#of gates: ', gateNum)
print('Calc # of gates: ', 4*n*l -2*l, 'gates')
#Convert the result from bin to int and compute the subtraction


#Print the results
#print('Expected result: ', int(max(np.dot(L,x)+b) - max(np.dot(M,x)+c)  )   )
#print('Garbled result:  ', res)
#print('Execution time:  ', time2-time1, 'seconds')
print('Size of Garbled circuit: ', str(sys.getsizeof(F)), 'bytes')


