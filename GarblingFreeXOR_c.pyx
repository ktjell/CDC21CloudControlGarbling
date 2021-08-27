# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:44:14 2021

@author: kst
"""

import numpy as np
import random
from hashlib import sha256
cimport numpy as cnp
ctypedef cnp.float64_t dtype_t
cimport cython
#random.seed(1)
cdef class garble:
    cdef int WIRE_BITLENGTH
    cdef int bitlen
    cdef int gate
    cdef long long R
    def __init__(self, int wire_bitlength, int bitlen):
        self.WIRE_BITLENGTH = wire_bitlength
        self.bitlen = bitlen
       
        
    def H(self,l,r,c=0,i=0):
        r'''
        Function that computes the hash of the sum of its inputs and returns the corresponding integer value.
        '''
        return int(sha256(str(l + r + c + i).encode('utf-8')).hexdigest(),16)
        
    def gg_MAX(self,tuple l, tuple r,tuple c,tuple k,int i):
        r'''
        Gate for choosing the maximum input based on the carry bit c. Like out = c*l + (1-c)*r . 
        '''
        cdef list out
        S0 = int(bin(k[0]) + self.WIRE_BITLENGTH*'0',2)     #The addition of the 128 zeros is what makes decryption possible.
        S1 = int(bin(k[1]) + self.WIRE_BITLENGTH*'0',2)        
       
        c1 = self.H(l[0], r[0], c[0], i)   ^ S0
        c2 = self.H(l[0], r[0], c[1], i)   ^ S0
        c3 = self.H(l[0], r[1], c[0], i)   ^ S0
        c4 = self.H(l[0], r[1], c[1], i)   ^ S1
        c5 = self.H(l[1], r[0], c[0], i)   ^ S1
        c6 = self.H(l[1], r[0], c[1], i)   ^ S0
        c7 = self.H(l[1], r[1], c[0], i)   ^ S1
        c8 = self.H(l[1], r[1], c[1], i)   ^ S1

        out = [c1,c2,c3,c4,c5,c6,c7,c8]
        random.shuffle(out)

 
        return out
    
    def MAXround(self,list l, list r):
        r'''
        Computes max{l,r}, where l and r are reprensented by bits (in fact l and r are keys of each bit of the input).
        The max is computed by first adding l with the 2's complement of r (using full adder gates) and the leading bit, R, of this computation
        determines the max. If R = 1 r is maximum and R=0, l is maximum.
        The output is the garbled gates and the keys of the output wires.
        '''
        #l and r are lists with bitlen number of keys.
        cdef int i
        cdef list gates = []
        cdef tuple and0 = (random.getrandbits(self.WIRE_BITLENGTH),)
        cdef tuple and1 = (np.bitwise_xor(and0[0], self.R),)
        gates.append(self.AND2(l[0], r[0][::-1], (and0[0], and1[0]), self.gate))
        self.gate += 1
        
        cdef list carry = [(and0[0],and1[0])]
        cdef list g
        cdef tuple s
        cdef tuple c
        for i in range(1,self.bitlen):
            g,s,c = self.FA0(l[i], r[i][::-1], carry[i-1])
            carry.append(c)
            gates.append(g)
        
        cdef tuple cc1 = self.XOR(carry[-1], carry[-2])
        cdef tuple cc2 = self.XOR(cc1, s)
       
        #Generate the gates for choosing the max of the left and right input
        cdef list out_keys0 = [random.getrandbits(self.WIRE_BITLENGTH) for i in range(self.bitlen)]
        cdef list out_keys1 = list(np.bitwise_xor(out_keys0, self.R))
        cdef list out_keys = [(k,j) for k,j in zip(out_keys0, out_keys1)]
        
        for i in range(self.bitlen):
            gates.append( self.gg_MAX(l[i], r[i], cc2, out_keys[i],  self.gate ) )
            self.gate+=1

        return gates, out_keys
    
    def AND(self,tuple l, tuple r, tuple o, int i):
        r'''
        Gate for choosing the maximum input based on the carry bit c. Like out = c*l + (1-c)*r . 
        '''
        cdef list out
        cdef tuple S0 = (int(bin(o[0]) + self.WIRE_BITLENGTH*'0',2) ,)    #The addition of the 128 zeros is what makes decryption possible.
        cdef tuple S1 = (int(bin(o[1]) + self.WIRE_BITLENGTH*'0',2) ,)        
       
        c1 = self.H(l[0], r[0], i)   ^ S0[0]
        c2 = self.H(l[0], r[1], i)   ^ S0[0]
        c3 = self.H(l[1], r[0], i)   ^ S0[0]
        c4 = self.H(l[1], r[1], i)   ^ S1[0]

        out = [c1,c2,c3,c4]
        random.shuffle(out)
        return out
    
    def AND2(self,tuple l, tuple r,tuple o,int i):
        r'''
        Gate for choosing the maximum input based on the carry bit c. Like out = c*l + (1-c)*r . 
        '''
        cdef list out
        cdef tuple S0 = (int(bin(o[0]) + self.WIRE_BITLENGTH*'0',2) ,)    #The addition of the 128 zeros is what makes decryption possible.
        cdef tuple S1 = (int(bin(o[1]) + self.WIRE_BITLENGTH*'0',2) ,)       
       
        c1 = self.H(l[0], r[0], i)   ^ S0[0]
        c2 = self.H(l[0], r[1], i)   ^ S1[0]
        c3 = self.H(l[1], r[0], i)   ^ S1[0]
        c4 = self.H(l[1], r[1], i)   ^ S1[0]

        out = [c1,c2,c3,c4]
        random.shuffle(out)
        return out
    

    def XOR(self,tuple l, tuple r):
        cdef tuple r0 = (np.bitwise_xor(l[0],r[0]),)
        return ( r0[0], np.bitwise_xor(r0[0],self.R) )

    def FA0(self,tuple l,tuple r , tuple c):
#        print(c)
        cdef tuple and0 = (random.getrandbits(self.WIRE_BITLENGTH),)
        cdef tuple and1 = (and0[0], np.bitwise_xor(and0[0], self.R))
        cdef tuple x1 = self.XOR(l,c)
#        print('x1', x1)
        cdef tuple x2 = self.XOR(r,c)
        cdef tuple x3 = self.XOR(x1,r)
        cdef list aGate = self.AND(x1, x2, and1, self.gate)
#        print('and0,and1', and0,and1)
        self.gate+=1
        cdef tuple x4 = self.XOR(and1,c)
        
        return aGate, x3, x4
        
    
    def ADDround(self, list l, list r):
        r'''
        Computes the garbled gates for adding two inputs. l and r are lists of (keys of) the bits of the inputs. Thus, this funtion calls a Full Adder
        gate for each pair of bits in the inputs. The output is the garbled gates and the output wires of the gate.
        '''
        cdef list gates = []
        cdef list sum_keys = []
        cdef tuple c = (0,0)
        cdef int i
        cdef list g
        cdef tuple s
#        print(l)
        for i in range(self.bitlen):
            g,s,c = self.FA0(l[i], r[i], c)
#            print('c',c)
            gates.append(g)
            sum_keys.append(s)

        return gates, sum_keys
        
    def garble(self, dict inputs):
        r'''
        Picks 2 random labels for each wire in the gate and "builds" the circuit up by using these labels and
        garbled versions of the gates.
        '''
        self.gate = 0
        cdef int n = len(inputs.keys())-1
        self.R = random.getrandbits(63)
        cdef int i,j        
        cdef list left_keys0 = [[random.getrandbits(self.WIRE_BITLENGTH) for i in range(self.bitlen)] for i in range(n)]
        cdef list left_keys = [[(ii, np.bitwise_xor(ii,self.R)) for ii in jj] for jj in left_keys0] 

        cdef list right_keys0 = [[random.getrandbits(self.WIRE_BITLENGTH) for i in range(self.bitlen)] for i in range(n)]
        cdef list right_keys = [[(ii,np.bitwise_xor(ii,self.R)) for ii in jj] for jj in right_keys0] 

        #left_keys are "flattened" because this is the format for (oblivious) transferring Alice's garbled input to Alice. 
        cdef list LEFTinp = []
        for ii in left_keys:
            LEFTinp.extend(ii)
        
        
        #Building up the garbled circuit by first adding the inputs pairwise
        cdef list gates = []              #Holds the gabled gates
        cdef list sum1_keys = []          #The output wires of the Full adder gates
        cdef list g,k
        for i in range(n):
            g,k = self.ADDround( left_keys[i], right_keys[i]) 
            gates.extend(g)
            sum1_keys.append(k)

        
        cdef list left  = sum1_keys[::2]
        cdef list right = sum1_keys[1::2]
        cdef list out_keys = []
        cdef list keys
        
        #Finding the MAX by pairwise comparison of the outputs of the addition. Assume n is power of 2
        #This is built like a tournament, so in first round n/2 values are compared and next n/4 and so on, so in total ceil(log(n)) rounds.
        
        for j in range(int(np.log2(n))):
            
            out_keys = []
            for i in range(len(right)):
                g, keys = self.MAXround(left[i],right[i])
                gates.append(g)
                out_keys.append(keys)
                
            left  = out_keys[::2]
            right = out_keys[1::2]   
        
        #### Adding the random number to mask the result from the evaluator
       
        right_k0 =[random.getrandbits(self.WIRE_BITLENGTH) for i in range(self.bitlen)] 
        right_k = [(ii,np.bitwise_xor(ii,self.R)) for ii in right_k0]
 
        g,k = self.ADDround(out_keys[0] , right_k) 
        gates.extend(g)
        out_keys = k
            
        right_keys.append(right_k)
#        print(len(right_keys))
#        The garbled version of Bob's input.
        cdef list Y = []          
        for i,ii in enumerate(right_keys):
            for j in range(self.bitlen):
                Y.append(ii[j][inputs[i][j]])
        

        return gates, Y, out_keys, LEFTinp #returns garbled curcuit, Bob's garbled input and decrypt info

#############################################################################################################
#############################################################################################################
        ## EVALUATION
        
    def Eval(self,list g,tuple l, tuple r, c=0,int i=0):
        r'''
        Takes as input a gate, the input wires to the gate and the gate number. It then evaluates the gate and decides 
        the "right" output by checking if t = 0. If t is never 0 something has gone wrong and the exection is aborted. 
        '''

        cdef tuple h = (self.H(l[0],r[0],c,i),)
        cdef tuple K    
        for i in range(len(g)):
            K = (g[i] ^ h[0],)
            if int(bin(K[0])[-self.WIRE_BITLENGTH:],2) == 0:
                return int(bin(K[0])[:-self.WIRE_BITLENGTH],2)
        
        raise ValueError('ABORT! Gate cannot be evaluated')

    def EvalADD(self, list g, list l, list r, list carry_keys = [0],int j=0):
        r'''
        g is a list of all garbled gates for computing the l+r. l and r are the keys for the bits in l and r. The function goes through
        the gates in g and evaluates them using the correct wires. output is the keys of the output wires.
        '''
        
        cdef list sum_keys = []
        cdef int i
        cdef tuple x1,x2
        for i in range(j, self.bitlen):
            x1 = (np.bitwise_xor(int(str(l[i])),carry_keys[i-j]),)
            x2 = (np.bitwise_xor(int(str(r[i])),carry_keys[i-j]),)
            sum_keys.append(np.bitwise_xor(x1[0],int(str(r[i]))))
            
            A1 = self.Eval(g[i], x1, x2, self.gate)
            self.gate+=1
            carry_keys.append(np.bitwise_xor(A1,carry_keys[i-j]))            
            
        
        return sum_keys, carry_keys

    
    def EvalMAX(self,list g, list l, list r):
        r'''
        g is a list of all garbled gates for computing the max{l,r}. l and r are the keys for the bits in l and r. The function goes through 
        the gates in g and evaluates them using the correct wires. output is the keys of the output wires.
        '''
        #l and r is lists with bitlen number of bits
        #gates is a list with the number of gates that is in the MAX-round (from Bob)
        
        #Full Adder part for adding left with minus right
                
        cdef tuple carry = (self.Eval(g[0],(l[0],),(r[0],),self.gate),)
        self.gate+=1
        cdef list s,c
        cdef tuple cc1, cc2
        cdef int i,j
        s,c = self.EvalADD(g[:self.bitlen], l, r, carry_keys = [carry[0]], j=1)
        cc1 = (np.bitwise_xor(c[-1], c[-2]),)
        cc2 = (np.bitwise_xor(cc1[0], s[-1]),)
        
        cdef list out_keys = []
        cdef list gates = g[self.bitlen:]
        for i in range(self.bitlen):
            out_keys.append(self.Eval(gates[i], (l[i],), (r[i],), cc2[0], self.gate ))
            self.gate+=1
        
        return out_keys
    

    
    def evaluate(self, list F, list Y, list X, list d, int n):
        r'''
        Takes the garbled circuit, the garbled input of the garbled, and the garbled input of the evaluator and decrypt information. Then the whole 
        garbled circuit is evaluated and finally the result is determined by using the decrypt info.
        If the result does not match with the decrypt information, something has gone wrong and the evaluator aborts the execution. 
        '''
        self.gate = 0
        #Evalaute the garbled circuit
        #FULL ADDER PART

        cdef list left  = []
        cdef list right = []
        cdef list gates = []
        cdef list out_keys
        cdef int i,j
        for i in range(n):
            left.append(X[(self.bitlen) * i : (self.bitlen) * (i+1) ])
            right.append(Y[(self.bitlen) * i : (self.bitlen) * (i+1) ])
            gates.append(F[(self.bitlen) * i : (self.bitlen) * (i+1) ])
#        print('LEFT',left)
        cdef list L2 = []
        for i in range(n):
            L2.append(self.EvalADD(gates[i], left[i], right[i], carry_keys = [0], j=0)[0])
            
        
        left = L2[::2]
        right= L2[1::2]
        cdef int GATE = self.bitlen*n
        cdef list keys
#        print('LEEEft', left)
        for j in range(int(np.log2(n))):
            
            out_keys = []
            for i in range(len(right)):
                keys = self.EvalMAX(F[GATE],left[i],right[i])
                GATE+=1
                out_keys.append(keys)
            
            left  = out_keys[::2]
            right = out_keys[1::2]
#        
        right = out_keys[0]
        left  = Y[-self.bitlen:]
        gates = F[-self.bitlen:]
        
        cdef list result = self.EvalADD(gates, left, right, carry_keys = [0], j=0)[0]
        
             
        #Decode the result
        cdef list res = []    
        for i in range(len(result)):
            if result[i] == d[i][0]:
                res.append(0)
            elif result[i] == d[i][1]:
                res.append(1)
            else:
                raise ValueError('ABORT, result could not be decoded.')
        
        
        cdef list res3 = res[::-1]
        
        return res3
















