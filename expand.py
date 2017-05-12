# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 21:09:24 2017

@author: liyi
"""

from sympy import *

def mean_sym(l):
    return sum(l)/len(l)

def pow2_sym(l):
    return [pow(i,2) for i in l]
            
def var_sym(l):
    M=mean_sym(l)
    L=[pow(i-M,2) for i in l]
    return mean_sym(L)            

r = Symbol('r')
s = Symbol('s')

r=0
s=1

f3 = Symbol('f3')
f4 = Symbol('f4')
f5 = Symbol('f5')
f6 = Symbol('f6')
f7 = Symbol('f7')
f88 = Symbol('f88')
f89 = Symbol('f89')
f90 = Symbol('f90')
f91 = Symbol('f91')
f92 = Symbol('f92')

x4=r*f3+f4+r*f5
x5=r*f4+f5+r*f6
x6=r*f5+f6+r*f7

x89=r*f88+f89+r*f90+s*r*f5
x90=r*f89+f90+r*f91+s*f5
x91=r*f90+f91+r*f92+s*r*f5

l1=[x5*x90,x4*x90,x6*x90,x5*x89,x5*x91,x4*x89,x4*x91,x6*x89,x6*x91]


var1=var_sym(l1)
    
simplify(expand(var1))

simplify(expand(mean_sym(l1)))




