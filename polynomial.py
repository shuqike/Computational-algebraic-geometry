from copy import deepcopy
from math import remainder
from pyclbr import Function
import numpy as np
from monomial import *

def cmp_to_key(order, monomials):
    class K:
        def __init__(self, obj, *args) -> None:
            self.obj = obj
        def __lt__(self, other):
            return order(monomials[self.obj], monomials[other.obj]) == 'greater'
        def __gt__(self, other):
            return order(monomials[self.obj], monomials[other.obj]) == 'less'
        def __eq__(self, other):
            return order(monomials[self.obj], monomials[other.obj]) == 'equal'
        def __le__(self, other):
            tmp = order(monomials[self.obj], monomials[other.obj])
            return tmp == 'greater' or tmp == 'equal'
        def __ge__(self, other):
            tmp = order(monomials[self.obj], monomials[other.obj])
            return tmp == 'less' or tmp == 'equal'
        def __ne__(self, other):
            return order(monomials[self.obj], monomials[other.obj]) != 'equal'
    return K

class Polynomial(object):
    def __init__(self, M, C) -> None:
        self.monomials = tuple(M)
        self.coeff = tuple(C)
    
    def reduce(self):
        self.isZero = True
        new_monomials = []
        new_coeff = []
        n = len(self.monomials)
        ind = np.zeros(n)
        for i in range(n):
            if ind[i] == 1:
                continue
            coeff = self.coeff[i]
            ind[i] = 1
            for j in range(i + 1, n):
                if self.monomials[i] == self.monomials[j]:
                    coeff += self.coeff[j]
                    ind[j] = 1
            if coeff != 0:
                self.isZero = False
                new_monomials.append(self.monomials[i])
                new_coeff.append(coeff)
        self.monomials = tuple(new_monomials)
        self.coeff = tuple(new_coeff)
    
    def sort(self, order):
        self.reduce()
        n = len(self.monomials)
        ind = sorted(np.arange(n), key=cmp_to_key(order, self.monomials))
        new_monomials = []
        new_coeff = []
        for i in range(n):
            new_monomials.append(self.monomials[ind[i]])
            new_coeff.append(self.coeff[ind[i]])
        self.monomials = tuple(new_monomials)
        self.coeff = tuple(new_coeff)
    
    def __str__(self) -> str:
        res = ''
        n = len(self.monomials)
        for i in range(n):
            res += '+(' + str(self.coeff[i]) + ')' + self.monomials[i].__str__()
        return res
    
    def __add__(self, other):
        new_monomials = list(self.monomials) + list(other.monomials)
        new_coeff = list(self.coeff) + list(other.coeff)
        p = Polynomial(tuple(new_monomials), tuple(new_coeff))
        p.reduce()
        return p
    
    def __sub__(self, other):
        new_monomials = list(self.monomials) + list(other.monomials)
        new_coeff = list(self.coeff) + list(np.asarray(other.coeff) * -1)
        p = Polynomial(tuple(new_monomials), tuple(new_coeff))
        p.reduce()
        return p
    
    def __mul__(self, other):
        n = len(self.monomials)
        m = len(other.monomials)
        new_monomials = []
        new_coeff = []
        for i in range(n):
            for j in range(m):
                new_monomials.append(self.monomials[i] * other.monomials[j])
                new_coeff.append(self.coeff[i] * other.coeff[j])
        p = Polynomial(tuple(new_monomials), tuple(new_coeff))
        p.reduce()
        return p
    
    def is_zero(self):
        self.reduce()
        return self.isZero
    
    def is_constant(self):
        if self.is_zero():
            return True
        for i in self.monomials:
            if (np.asarray(i.power) == 0).all():
                return True
        return False
    
    def __eq__(self, other) -> bool:
        if self.is_zero() and other == 0:
            return True
        if type(other) != Polynomial:
            return False
        return self.monomials == other.monomials and self.coeff == other.coeff

def S_poly(p1:Polynomial, p2:Polynomial, order):
    p1.sort(order)
    p2.sort(order)
    x_gamma = mono_termwise_max(p1.monomials[0], p2.monomials[0])
    q1 = Polynomial([x_gamma / p1.monomials[0]], [1 / p1.coeff[0]])
    q2 = Polynomial([x_gamma / p2.monomials[0]], [1 / p2.coeff[0]])
    return q1 * p1 - q2 * p2

def division_alg(p:Polynomial, F, order):
    s = len(F)
    Q = [0 for _ in range(s)]
    remainder = 0
    while p.is_zero() == False:
        i = 0
        division_occurred = False
        while i < s and division_occurred == False:
            F[i].sort(order)
            p.sort(order)
            if mono_divisible(F[i].monomials[0], p.monomials[0]):
                tmp = Polynomial([p.monomials[0] / F[i].monomials[0]], [p.coeff[0] / F[i].coeff[0]])
                if Q[i] == 0:
                    Q[i] = tmp
                else:
                    Q[i] = Q[i] + tmp
                p = p - tmp * F[i]
                division_occurred = True
            else:
                i = i + 1
        if division_occurred == False:
            tmp = Polynomial([p.monomials[0]], [p.coeff[0]])
            if remainder == 0:
                remainder = tmp
            else:
                remainder = remainder + tmp
            p = p - tmp
    return Q, remainder

def poly_equiv(p1, p2):
    if type(p1) != Polynomial and type(p2) != Polynomial:
        if p1 != 0 and p2 != 0:
            return True
        if p1 == 0 and p2 == 0:
            return True
        return False
    if type(p1) != Polynomial:
        if p2.is_constant() == False:
            return False
        if p1 != 0 and p2.is_zero() == False:
            return True
        if p1 == 0 and p2.is_zero():
            return True
        return False
    if type(p2) != Polynomial:
        if p1.is_constant() == False:
            return False
        if p2 != 0 and p1.is_zero() == False:
            return True
        if p2 == 0 and p1.is_zero():
            return True
        return False
    if p1.is_zero() != p2.is_zero():
        return False
    factor = 0
    n = len(p1.monomials)
    m = len(p2.monomials)
    for j in range(m):
        if p2.monomials[j] == p1.monomials[0]:
            factor = p1.coeff[0] / p2.coeff[j]
    if factor == 0:
        return False
    for i in range(n):
        match = False
        for j in range(m):
            if p2.monomials[j] == p1.monomials[i]:
                if p2.coeff[j] * factor != p1.coeff[i]:
                    return False
                match = True
                break
        if match == False:
            return False
    return True

def Buchberger_alg(F, order):
    s = len(F)
    G = list(deepcopy(F))
    cnt = 0
    while True:
        cnt += 1
        newG = deepcopy(G)
        length = len(newG)
        expanded = False
        for i in range(length):
            for j in range(i + 1, length):
                Q, r = division_alg(S_poly(newG[i], newG[j], order), newG, order)
                if r != 0 or (type(r) == Polynomial and r.is_zero() == False):
                    equiv_exist = False
                    for g in G:
                        if poly_equiv(g, r):
                            equiv_exist = True
                            break
                    if equiv_exist == False:
                        expanded = True
                        G.append(r)
        if expanded == False:
            break
    return G

def test1():
    '''
    Test reduce
    '''
    m1 = Monomial((1,2,3))
    c1 = 1
    m2 = Monomial((1,2,3))
    c2 = 2
    m3 = Monomial((0,0,1))
    c3 = 5
    Ms = [m1, m2, m3]
    Cs = [c1, c2, c3]
    p = Polynomial(Ms, Cs)
    p.reduce()
    print(p)

def test2():
    '''
    Test sort
    Correct output:
    +(-5)[3, 0, 0]+(7)[2, 0, 2]+(4)[1, 2, 1]+(4)[0, 0, 2]
    +(7)[2, 0, 2]+(4)[1, 2, 1]+(-5)[3, 0, 0]+(4)[0, 0, 2]
    +(4)[1, 2, 1]+(7)[2, 0, 2]+(-5)[3, 0, 0]+(4)[0, 0, 2]
    '''
    m1 = Monomial((3,0,0))
    m2 = Monomial((2,0,2))
    m3 = Monomial((1,2,1))
    m4 = Monomial((0,0,2))
    Ms = [m1, m2, m3, m4]
    Cs = [-5, 7, 4, 4]
    p = Polynomial(Ms, Cs)
    p.sort(lex_order)
    print(p)
    p.sort(grlex_order)
    print(p)
    p.sort(grevlex_order)
    print(p)

def test3():
    '''
    Test S_poly
    Correct output:
    +(-1.0)[3, 3]+(1.0)[2, 0]+(-0.3333333333333333)[0, 3]
    '''
    m1 = Monomial((3,2))
    m2 = Monomial((2,3))
    m3 = Monomial((1,0))
    Ms = [m1, m2, m3]
    Cs = [1, -1, 1]
    p1 = Polynomial(Ms, Cs)
    m1 = Monomial((4,1))
    m2 = Monomial((0,2))
    Ms = [m1, m2]
    Cs = [3, 1]
    p2 = Polynomial(Ms, Cs)
    print(S_poly(p1, p2, grlex_order))

def test4():
    '''
    Test division algorithm
    Correct output:
    +(1.0)[1, 0]+(1.0)[0, 0]
    +(1.0)[1, 0]
    +(2.0)[1, 0]+(1.0)[0, 0]
    '''
    m1 = Monomial((0,2))
    m2 = Monomial((0,0))
    Ms = [m1, m2]
    Cs = [1, -1]
    f1 = Polynomial(Ms, Cs)
    m1 = Monomial((1,1))
    m2 = Monomial((0,0))
    Ms = [m1, m2]
    Cs = [1, -1]
    f2 = Polynomial(Ms, Cs)
    F = (f1, f2)
    m1 = Monomial((2,1))
    m2 = Monomial((1,2))
    m3 = Monomial((0,2))
    Ms = [m1, m2, m3]
    Cs = [1, 1, 1]
    p = Polynomial(Ms, Cs)
    Q, remainder = division_alg(p, F, lex_order)
    print(Q[0])
    print(Q[1])
    print(remainder)

def test5():
    '''
    Test Buchberger’s Algorithm
    Correct output:
    '''
    m1 = Monomial((3,0))
    m2 = Monomial((1,1))
    Ms = [m1, m2]
    Cs = [1, -2]
    f1 = Polynomial(Ms, Cs)
    m1 = Monomial((2,1))
    m2 = Monomial((0,2))
    m3 = Monomial((1,0))
    Ms = [m1, m2, m3]
    Cs = [1, -2, 1]
    f2 = Polynomial(Ms, Cs)
    F = [f1, f2]
    G = Buchberger_alg(F, grlex_order)
    for g in G:
        print(g)
