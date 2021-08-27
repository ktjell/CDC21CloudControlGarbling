# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:14:29 2021

@author: kst
"""
import numpy as np

class ElGamal:
    def __init__(self, G, g, q):
        self.G = G
        self.g = g
        self.q = q
        
    def Gen(self, sk):
        r'''
        Public key generation algortihm of the ElGamal crypto scheme.
        '''
        h = self.G( pow( int(str(self.g)), int(sk), self.G.p))                      
        return self.g, h
    
    def OGen(self):
        r'''
        Oblivious public key generation, which works for the group Z_p, where p = 2q + 1, where q and p are primes.
        '''
        r = np.random.randint(0,2**31)  #Just chosing a random number from Z_p.
        return (self.g, self.G(pow(r,2,self.G.p)))
    
    def Enc(self, m, pk):
        r'''
        Encryption algortihm of the ElGamal crypto scheme.
        '''
        r = np.random.randint(0, 2**31 )
        
        return (self.G(pow( int(str( pk[0] )) ,r, self.G.p)), self.G(m*pow( int(str(pk[1])) ,r, self.G.p)))
    
    def Dec(self, c, sk):
        r'''
        Decryption algortihm of the ElGamal crypto scheme.
        '''

        Cinv = int(str(c[0].inverse()))
        m =c[1] * pow(Cinv, sk, self.G.p)

        return m 
