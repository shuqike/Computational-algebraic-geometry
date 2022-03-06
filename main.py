from sympy import Poly
from monomial import Monomial
from polynomial import *

'''a'''
def problem_a(order):
    m1 = Monomial((2,1))
    m2 = Monomial((0,0))
    Ms = [m1, m2]
    Cs = [1, -1]
    f1 = Polynomial(Ms, Cs)
    m1 = Monomial((1,2))
    m2 = Monomial((1,0))
    Ms = [m1, m2]
    Cs = [1, -1]
    f2 = Polynomial(Ms, Cs)
    F = [f1, f2]
    G = Buchberger_alg(F, order)
    for g in G:
        print(g)

print('--Problem a--')
print('lex order:')
problem_a(lex_order)
print('grlex order:')
problem_a(grlex_order)
print('--end')

'''b'''
def problem_b(order):
    m1 = Monomial((2,0))
    m2 = Monomial((0,1))
    Ms = [m1, m2]
    Cs = [1, 1]
    f1 = Polynomial(Ms, Cs)
    m1 = Monomial((4,0))
    m2 = Monomial((2,1))
    m3 = Monomial((0,2))
    m4 = Monomial((0,0))
    Ms = [m1, m2, m3, m4]
    Cs = [1, 2, 1, 3]
    f2 = Polynomial(Ms, Cs)
    F = [f1, f2]
    G = Buchberger_alg(F, order)
    for g in G:
        print(g)

print('--Problem b--')
print('lex order:')
problem_b(lex_order)
print('grlex order:')
problem_b(grlex_order)
print('--end')

'''c'''
def problem_c(order):
    m1 = Monomial((1,0,0))
    m2 = Monomial((0,0,4))
    Ms = [m1, m2]
    Cs = [1, -1]
    f1 = Polynomial(Ms, Cs)
    m1 = Monomial((0,1,0))
    m2 = Monomial((0,0,5))
    Ms = [m1, m2]
    Cs = [1, -1]
    f2 = Polynomial(Ms, Cs)
    F = [f1, f2]
    G = Buchberger_alg(F, order)
    for g in G:
        print(g)

print('--Problem c--')
print('lex order:')
problem_c(lex_order)
print('grlex order:')
problem_c(grlex_order)
print('--end')
