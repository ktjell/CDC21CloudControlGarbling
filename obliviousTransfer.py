# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:42:09 2021

@author: kst
"""
import numpy as np

class OT:
    def __init__(self, bitlen, elg, mod):
        self.bitlen = bitlen
        self.elg = elg
        self.mod = mod
        
    def choose(self, inputs):
        r'''
        First fase of OT. Alice choose 3 real public keys and 3 fake keys and sends them to Bob.
        '''
        
        w1 = self.bitlen*len(inputs.keys())
        self.sk = [np.random.randint(1,self.mod, dtype = np.int64) for i in range(w1)]      #The "true" secret keys choosen by Alice
       
        pkb = [self.elg.Gen(i) for i in self.sk]               #The true public key, generated with the true secret key.
        self.pk = []                      
        for i in range(w1):
            self.pk.append(self.elg.OGen())                    #Creating the falsh public keys.
        
        
        indices = list(inputs.values())
        indices = np.array(indices).flatten()
        for i in range(w1):
            self.pk.insert(indices[i]+2*i, pkb[i])      #Insert the true public key at the index corresponding to Alice's input.
        
        return self.pk   
        
    def transfer(self, pk, mes):
        r'''
        The middle fase of OT. Bob encrypts the labels for the wires of Alice's inputs and sends to Alice.
        '''
        c = []
        for i in range(len(mes)):
            c.append((self.elg.Enc(mes[i][0], pk[2*i]), self.elg.Enc(mes[i][1], pk[2*i+1])))
        return c
    

    def retrieve(self, c, inputs):
        r'''
        Alice decrypts the message from Bob and gets the garbled version of her input to the circuit.
        '''
        w1 = self.bitlen*len(inputs.keys())
        indices = list(inputs.values())
        indices = np.array(indices).flatten()
        mes = []
        for i in range(w1):
            mes.append( self.elg.Dec(c[i][indices[i]], int(self.sk[i])) )
        
        return mes
    