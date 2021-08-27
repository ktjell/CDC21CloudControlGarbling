# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 16:44:14 2021

@author: kst
"""

import numpy as np
import random
from hashlib import sha256

class garble:
    def __init__(self, wire_bitlength, bitlen):
        self.WIRE_BITLENGTH = wire_bitlength
        self.bitlen = bitlen
       
        
    def H(self,l,r,c=0,i=0):
        r'''
        Function that computes the hash of the sum of its inputs and returns the corresponding integer value.
        '''
        return int(sha256(str(l + r + c + i).encode('utf-8')).hexdigest(),16)
        
    def gg_MAX(self,l, r, c, k, i):
        r'''
        Gate for choosing the maximum input based on the carry bit c. Like out = c*l + (1-c)*r . 
        '''
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

    def gg_OVERflow(self,l, r, c, k, i):
        r'''
        Sometimes "overflow" occurs when adding input 1 with - input 2 for computing max. This gate corrects that carry bit.
        '''
        S0 = int(bin(k[0]) + self.WIRE_BITLENGTH*'0',2)     #The addition of the 128 zeros is what makes decryption possible.
        S1 = int(bin(k[1]) + self.WIRE_BITLENGTH*'0',2)        
       
        c1 = self.H(l[0], r[0], c[0], i)   ^ S0
        c2 = self.H(l[0], r[0], c[1], i)   ^ S1
        c3 = self.H(l[0], r[1], c[0], i)   ^ S1
        c4 = self.H(l[0], r[1], c[1], i)   ^ S0
        c5 = self.H(l[1], r[0], c[0], i)   ^ S1
        c6 = self.H(l[1], r[0], c[1], i)   ^ S0
        c7 = self.H(l[1], r[1], c[0], i)   ^ S0
        c8 = self.H(l[1], r[1], c[1], i)   ^ S1

        out = [c1,c2,c3,c4,c5,c6,c7,c8]
        random.shuffle(out)

 
        return out

    def gg_FA2(self,l,r,c,k, K, i):
        r'''
        Computes the addition out = inp1 + (-inp2). The first half is returned for the half adder (when the carry is known to be 0)
        '''
        S0 = int(bin(k[0]) + self.WIRE_BITLENGTH*'0',2)
        S1 = int(bin(k[1]) + self.WIRE_BITLENGTH*'0',2)
        
        C0 = int(bin(K[0]) + self.WIRE_BITLENGTH*'0',2)
        C1 = int(bin(K[1]) + self.WIRE_BITLENGTH*'0',2)
        
        
        c1 = self.H(l[0], r[0], c[1], i)   ^ S0
        d1 = self.H(l[0], r[0], c[1], i+1) ^ C1
        c2 = self.H(l[0], r[1], c[1], i)   ^ S1
        d2 = self.H(l[0], r[1], c[1], i+1) ^ C0
        c3 = self.H(l[1], r[0], c[1], i)   ^ S1
        d3 = self.H(l[1], r[0], c[1], i+1) ^ C1       
        c4 = self.H(l[1], r[1], c[1], i)   ^ S0
        d4 = self.H(l[1], r[1], c[1], i+1) ^ C1       
        
        if c[0] == None:         #For the first bits to be added, the carry bit will be one, therefor we dont need
            c = [c1,c2,c3,c4]   #to compute the result for all c that is zero.
            d = [d1,d2,d3,d4]
            random.shuffle(c)
            random.shuffle(d)
            return c, d
 
        c5 = self.H(l[0], r[0], c[0], i)   ^ S1
        d5 = self.H(l[0], r[0], c[0], i+1) ^ C0    
        c6 = self.H(l[0], r[1], c[0], i)   ^ S0
        d6 = self.H(l[0], r[1], c[0], i+1) ^ C0
        c7 = self.H(l[1], r[0], c[0], i)   ^ S0
        d7 = self.H(l[1], r[0], c[0], i+1) ^ C1
        c8 = self.H(l[1], r[1], c[0], i)   ^ S1
        d8 = self.H(l[1], r[1], c[0], i+1) ^ C0
                
        c = [c1,c2,c3,c4,c5,c6,c7,c8]
        d = [d1,d2,d3,d4,d5,d6,d7,d8]
        
        random.shuffle(c)
        random.shuffle(d)
        return c, d

        
    def gg_FA(self,l,r,c,k, K, i):
        r'''
        Computes a normal addition. out = inp1 + inp2. The first half is returned for the half adder (when the carry is known to be 0)
        '''
        S0 = int(bin(k[0]) + self.WIRE_BITLENGTH*'0',2)
        S1 = int(bin(k[1]) + self.WIRE_BITLENGTH*'0',2)
        
        C0 = int(bin(K[0]) + self.WIRE_BITLENGTH*'0',2)
        C1 = int(bin(K[1]) + self.WIRE_BITLENGTH*'0',2)
        
        
        c1 = self.H(l[0], r[0], c[0], i)   ^ S0
        d1 = self.H(l[0], r[0], c[0], i+1) ^ C0
        c2 = self.H(l[0], r[1], c[0], i)   ^ S1
        d2 = self.H(l[0], r[1], c[0], i+1) ^ C0
        c3 = self.H(l[1], r[0], c[0], i)   ^ S1
        d3 = self.H(l[1], r[0], c[0], i+1) ^ C0       
        c4 = self.H(l[1], r[1], c[0], i)   ^ S0
        d4 = self.H(l[1], r[1], c[0], i+1) ^ C1       
        
        if len(c) == 1:         #For the first bits to be added, the carry bit will be zero, there for we dont need
            c = [c1,c2,c3,c4]   #to compute the result for all c that is 1.
            d = [d1,d2,d3,d4]
            random.shuffle(c)
            random.shuffle(d)
            return c, d
 
        c5 = self.H(l[0], r[0], c[1], i)   ^ S1
        d5 = self.H(l[0], r[0], c[1], i+1) ^ C0    
        c6 = self.H(l[0], r[1], c[1], i)   ^ S0
        d6 = self.H(l[0], r[1], c[1], i+1) ^ C1
        c7 = self.H(l[1], r[0], c[1], i)   ^ S0
        d7 = self.H(l[1], r[0], c[1], i+1) ^ C1
        c8 = self.H(l[1], r[1], c[1], i)   ^ S1
        d8 = self.H(l[1], r[1], c[1], i+1) ^ C1
                
        c = [c1,c2,c3,c4,c5,c6,c7,c8]
        d = [d1,d2,d3,d4,d5,d6,d7,d8]
        
        random.shuffle(c)
        random.shuffle(d)
        return c, d


    
    def MAXround(self,l,r):
        r'''
        Computes max{l,r}, where l and r are reprensented by bits (in fact l and r are keys of each bit of the input).
        The max is computed by first adding l with the 2's complement of r (using full adder gates) and the leading bit, R, of this computation
        determines the max. If R = 1 r is maximum and R=0, l is maximum.
        The output is the garbled gates and the keys of the output wires.
        '''
        #l and r are lists with bitlen number of keys.
        
        gates = []

        #Generating Full adder gates for adding the left input with the complement of the right input.
        
        carry2_keys = [(random.getrandbits(128), random.getrandbits(128)) for i in range(self.bitlen)]
        sum2_keys   = [(random.getrandbits(128), random.getrandbits(128)) for i in range(self.bitlen)]
        
        gates.extend(self.gg_FA2(l[0], r[0], [None,0], sum2_keys[0], carry2_keys[0], self.gate))
        self.gate+=2
        
        for i in range(1,self.bitlen):
            gates.extend( self.gg_FA2(l[i], r[i], carry2_keys[i-1], sum2_keys[i], carry2_keys[i], self.gate) )
            self.gate+=2
        
          
#        Generating gates for compensating for overflow in the previous full adder. 
        R_keys   = (random.getrandbits(128), random.getrandbits(128))
        gates.append(self.gg_OVERflow(carry2_keys[-2], carry2_keys[-1], sum2_keys[-1], R_keys, self.gate))
        self.gate+=1
#        R_keys = sum2_keys[-1]
                
        
        #Generate the gates for choosing the max of the left and right input
        out_keys = [(random.getrandbits(128), random.getrandbits(128)) for i in range(self.bitlen)]
        for i in range(self.bitlen):
            gates.append( self.gg_MAX(l[i], r[i], R_keys, out_keys[i],  self.gate ) )
            self.gate+=1

        return gates, out_keys

    
    def ADDround(self,l,r):
        r'''
        Computes the garbled gates for adding two inputs. l and r are lists of (keys of) the bits of the inputs. Thus, this funtion calls a Full Adder
        gate for each pair of bits in the inputs. The output is the garbled gates and the output wires of the gate.
        '''
        sum_keys = [(random.getrandbits(self.WIRE_BITLENGTH), random.getrandbits(self.WIRE_BITLENGTH)) for i in range(self.bitlen)]
        carry_keys = [(random.getrandbits(self.WIRE_BITLENGTH), random.getrandbits(self.WIRE_BITLENGTH)) for i in range(self.bitlen)]
        gates = []
        gates.extend(self.gg_FA(l[0],r[0],[0],sum_keys[0], carry_keys[0], self.gate))
        self.gate+=2
        for i in range(1,self.bitlen):
            gates.extend(self.gg_FA(l[i], r[i], carry_keys[i-1], sum_keys[i], carry_keys[i], self.gate ))
            self.gate+=2
        
#        gates.append(self.gg_corCarry( l[-1], r[-1], carry_keys[-1], sum_keys[-1], self.gate))
#        self.gate += 1

        return gates, sum_keys
        
    def garble(self, inputs):
        r'''
        Picks 2 random labels for each wire in the gate and "builds" the circuit up by using these labels and
        garbled versions of the gates.
        '''
        self.gate = 0
        num_bits = self.bitlen
        n = len(inputs.keys())
                
        #Generate the keys for the inputs. left_keys are Alice input and right:keys are Bob' input.
        left_keys  = [[(random.getrandbits(self.WIRE_BITLENGTH), random.getrandbits(self.WIRE_BITLENGTH)) for i in range(num_bits)] for i in range(n)]
        right_keys = [[(random.getrandbits(self.WIRE_BITLENGTH), random.getrandbits(self.WIRE_BITLENGTH)) for i in range(num_bits)] for i in range(n)]
        
        #left_keys are "flattened" because this is the format for (oblivious) transferring Alice's garbled input to Alice. 
        LEFTinp=[]
        for i in left_keys:
            LEFTinp.extend(i)
        
        #Building up the garbled circuit by first adding the inputs pairwise
        gates = []              #Holds the gabled gates
        sum1_keys = []          #The output wires of the Full adder gates
        for i in range(n):
            g,k = self.ADDround( left_keys[i], right_keys[i]) 
            gates.extend(g)
            sum1_keys.append(k)

        
        left  = sum1_keys[::2]
        right = sum1_keys[1::2]
        
        
        #Finding the MAX by pairwise comparison of the outputs of the addition.
        #This is built like a tournament, so in first round n/2 values are compared and next n/4 and so on, so in total ceil(log(n)) rounds.
        # The flag is used to keep track of the excess value if n is uneven. In this case, the extra value is compared in the next round
        # in the block "if flag:".
        
        flag = 0
        carry = 0
        
        for j in range(int(np.ceil(np.log2(n)))):
           
            out_keys = []
            
            for i in range(len(right)):
                g, keys = self.MAXround(left[i],right[i])
                gates.append(g)
                out_keys.append(keys)
            
            if flag:    #compare now the excess value from the previous round to the last left element of this round.
                g,keys = self.MAXround(left[-1],carry)
                gates.append(g)
                out_keys.append(keys)
                
            flag = n%2
            n = int(n/2) 
            if flag:    
                carry = left[-1] #keep the excess value that has not yet been compared
            
            left  = out_keys[::2]
            right = out_keys[1::2]   
            

        #The garbled version of Bob's input.
        Y = []          
        for i,ii in enumerate(right_keys):
            for j in range(num_bits):
                Y.append(ii[j][inputs[i][j]])
        

        return gates, Y, out_keys[0], LEFTinp #returns garbled curcuit, Bob's garbled input and decrypt info

#############################################################################################################
#############################################################################################################
        ## EVALUATION
        
    def Eval(self,g,l,r,c=0,i=0):
        r'''
        Takes as input a gate, the input wires to the gate and the gate number. It then evaluates the gate and decides 
        the "right" output by checking if t = 0. If t is never 0 something has gone wrong and the exection is aborted. 
        '''

        h = self.H(l,r,c,i)
            
        for i in range(len(g)):
            K = g[i] ^ h
            if int(bin(K)[-self.WIRE_BITLENGTH:],2) == 0:
                return int(bin(K)[:-self.WIRE_BITLENGTH],2)
        
        raise ValueError('ABORT! Gate cannot be evaluated')
        
    
    def EvalMAX(self,g,l,r):
        r'''
        g is a list of all garbled gates for computing the max{l,r}. l and r are the keys for the bits in l and r. The function goes through 
        the gates in g and evaluates them using the correct wires. output is the keys of the output wires.
        '''
        #l and r is lists with bitlen number of bits
        #gates is a list with the number of gates that is in the MAX-round (from Bob)
        
        #Full Adder part for adding left with minus right
        gN = 0#2*self.bitlen
        sum2_keys   = [self.Eval(g[gN],   l[0], r[0],0, self.gate)]
        carry2_keys = [self.Eval(g[gN+1], l[0], r[0],0, self.gate+1) ]
        self.gate += 2
        gN += 2        
        for i in range(1,self.bitlen):
            sum2_keys.append(  self.Eval(g[gN],  l[i], r[i], carry2_keys[i-1], self.gate ))
            carry2_keys.append(self.Eval(g[gN+1],l[i], r[i], carry2_keys[i-1], self.gate+1) )
            gN += 2
            self.gate += 2
            
        #Compensating for overflow
        R_keys = self.Eval(g[gN], carry2_keys[-2], carry2_keys[-1], sum2_keys[-1], self.gate)
        gN+=1
        self.gate += 1
#        R_keys = sum2_keys[-1]
        #Choosing MAX of left and right based on the "sign of the addition".
        out_keys = []
        for i in range(self.bitlen):
            out_keys.append(self.Eval(g[gN], l[i], r[i], R_keys, self.gate ))
            gN+=1
            self.gate+=1
        
        return out_keys
    
    def EvalADD(self, g,l,r):
        r'''
        g is a list of all garbled gates for computing the l+r. l and r are the keys for the bits in l and r. The function goes through
        the gates in g and evaluates them using the correct wires. output is the keys of the output wires.
        '''
        sum_keys   = [self.Eval(g[0], l[0], r[0], 0, self.gate)]
        carry_keys = [self.Eval(g[1], l[0], r[0], 0, self.gate+1) ]
        self.gate += 2
        gN = 2
        for i in range(1,self.bitlen):
            sum_keys.append(  self.Eval(g[gN],   l[i], r[i], carry_keys[i-1], self.gate ) )
            carry_keys.append(self.Eval(g[gN+1], l[i], r[i], carry_keys[i-1], self.gate+1 ) )
            self.gate += 2
            gN += 2
        
#        sum_keys.append(self.Eval(g[-1], l[-1], r[-1], carry_keys[-1], self.gate ) )
#        self.gate +=1
        return sum_keys
        
    
    def evaluate(self, F, Y, X, d, n):
        r'''
        Takes the garbled circuit, the garbled input of the garbled, and the garbled input of the evaluator and decrypt information. Then the whole 
        garbled circuit is evaluated and finally the result is determined by using the decrypt info.
        If the result does not match with the decrypt information, something has gone wrong and the evaluator aborts the execution. 
        '''
        self.gate = 0
        #Evalaute the garbled circuit
        #FULL ADDER PART

        left = []
        right = []
        gates = []
        for i in range(n):
            left.append(X[(self.bitlen) * i : (self.bitlen) * (i+1) ])
            right.append(Y[(self.bitlen) * i : (self.bitlen) * (i+1) ])
            gates.append(F[(2*self.bitlen) * i : (2*self.bitlen) * (i+1) ])

        L2 = []
        for i in range(n):
            L2.append(self.EvalADD(gates[i], left[i], right[i]))
            
        
        left = L2[::2]
        right= L2[1::2]
        
        #MAX part
        flag = 0
        GATE = self.gate
        carry = 0
        for j in range(int(np.ceil(np.log2(n)))):
            out_keys = []
            for i in range(len(right)):
                keys = self.EvalMAX(F[GATE],left[i],right[i])
                GATE+=1
                out_keys.append(keys)
            
            if flag:
                keys = self.EvalMAX(F[GATE], left[-1],carry)
                GATE+=1
                out_keys.append(keys)
                
            flag = n%2
            n = int(n/2)
            if flag:
                carry = left[-1]
            
            left  = out_keys[::2]
            right = out_keys[1::2]
        
        result = left[0]
             
        #Decode the result
        res = []    
        for i in range(len(result)):
            if result[i] == d[i][0]:
                res.append(0)
            elif result[i] == d[i][1]:
                res.append(1)
            else:
                raise ValueError('ABORT, result could not be decoded.')
        
        
        res3 = res[::-1]
        
        return res3
















