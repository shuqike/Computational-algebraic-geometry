import numpy as np

class Monomial(object):
    def __init__(self, power) -> None:
        self.power = np.asarray(power)
    
    def total_degree(self):
        return sum(self.power)
    
    def __mul__(self, other):
        return Monomial(self.power + other.power)
    
    def __truediv__(self, other):
        return Monomial(self.power - other.power)
    
    def __str__(self) -> str:
        x = list(self.power)
        return x.__str__()
    
    def __eq__(self, other):
        return (self.power == other.power).all()

def lex_order(m1:Monomial, m2:Monomial):
    x = m1.power - m2.power
    for i in x:
        if i < 0:
            return 'less'
        if i > 0:
            return 'greater'
    return 'equal'

def invlex_order(m1:Monomial, m2:Monomial):
    x = m1.power - m2.power
    n = len(x)
    for i in range(n):
        if x[n-i-1] < 0:
            return 'less'
        if x[n-i-1] > 0:
            return 'greater'
    return 'equal'

def grlex_order(m1:Monomial, m2:Monomial):
    t1 = m1.total_degree()
    t2 = m2.total_degree()
    if t1 < t2:
        return 'less'
    if t1 > t2:
        return 'greater'
    return lex_order(m1, m2)

def grevlex_order(m1:Monomial, m2:Monomial):
    t1 = m1.total_degree()
    t2 = m2.total_degree()
    if t1 < t2:
        return 'less'
    if t1 > t2:
        return 'greater'
    return invlex_order(m2, m1)

def mono_termwise_max(m1:Monomial, m2:Monomial):
    power = []
    if len(m1.power) != len(m2.power):
        print('||---monomial termwise max failed because of non-equal length---||')
        return None
    for i, j in zip(m1.power, m2.power):
        power.append(max(i, j))
    return Monomial(power)

def mono_divisible(m1:Monomial, m2:Monomial):
    power = m2.power - m1.power
    return (power >= 0).all()

def mono_test1():
    ''' All output should be 'greater' '''
    print(lex_order(Monomial((1,2,0)), Monomial((0,3,4))))
    print(lex_order(Monomial((3,2,4)), Monomial((3,2,1))))
    print(grlex_order(Monomial((1,2,3)), Monomial((3,2,0))))
    print(grlex_order(Monomial((1,2,4)), Monomial((1,1,5))))
    print(grevlex_order(Monomial((4,7,1)), Monomial((4,2,3))))
    print(grevlex_order(Monomial((1,5,2)), Monomial((4,1,3))))

def mono_test2():
    '''
    Correct output:
    True
    False
    '''
    print( mono_divisible( Monomial((1,1,1)), Monomial((1,2,3)) ) )
    print( mono_divisible( Monomial((1,3,1)), Monomial((1,2,3)) ) )
