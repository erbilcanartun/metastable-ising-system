import numpy as np
from scipy.misc import derivative

class RenormalizationGroup:

    def __init__(self):

        self.b = 2 # Decimation
        self.d = 3 # Dimension
        
        self.m = self.b**(self.d - 1) # Bond-moving multiplier
        
        self.eigen = self.b ** self.d
        self.neigen = self.b ** (-self.d)

        self.Jc = 0.065266294981 # Critical point

    def _pderivative(self, function, variable=0, point=[]):
        arguments = point[:]
        def wraps(x):
            arguments[variable] = x
            return function(*arguments)
        return derivative(wraps, point[variable], dx=np.sqrt(np.finfo(float).eps))

    # Recursion relations (..._1) for Hamiltonian : J.(SiSj) + H.(Si)
    def _lnR1_1(self, j, h):

        e1 = 2*j + h
        e2 = -2*j - h
        emax = np.amax([e1, e2])

        return emax + np.log(np.exp(e1 - emax) + np.exp(e2 - emax))

    def _lnR2_1(self, j, h):

        e1 = -2*j + h
        e2 = 2*j - h
        emax = np.amax([e1, e2])

        return emax + np.log(np.exp(e1 - emax) + np.exp(e2 - emax))

    def _lnR3_1(self, j, h):

        e1 = h
        e2 = -h
        emax = np.amax([e1, e2])

        return emax + np.log(np.exp(e1 - emax) + np.exp(e2 - emax))

    def J_1(self, interaction, field):
        j = self.m * interaction
        h = field
        return (1/4) * (self._lnR1_1(j, h) + self._lnR2_1(j, h) - 2*self._lnR3_1(j, h))

    def H_1(self, interaction, field):
        j = self.m * interaction
        h = field
        return (1/4) * (self._lnR1_1(j, h) - self._lnR2_1(j, h))

    def G_1(self, interaction, field):
        j = self.m * interaction
        h = field
        return (1/4) * (self._lnR1_1(j, h) + self._lnR2_1(j, h) + 2*self._lnR3_1(j, h))

    def _recursion_matrix_1(self, j, h):

        X = self._pderivative(self.G_1, 0, [j, h])
        Y = self._pderivative(self.J_1, 0, [j, h])
        Z = self._pderivative(self.H_1, 0, [j, h])

        A = self._pderivative(self.G_1, 1, [j, h])
        B = self._pderivative(self.J_1, 1, [j, h])
        C = self._pderivative(self.H_1, 1, [j, h])

        return np.array([[self.eigen, X,  self.eigen * A],
                         [0,          Y,  self.eigen * B],
                         [0,          Z,  self.eigen * C]])

    # Recursion relations for Hamiltonian : J.(SiSj) + H.(Si+Sj)
    def _lnR1(self, j, h):

        e1 = 2*j + 4*h
        e2 = -2*j
        emax = np.amax([e1, e2])

        return emax + np.log(np.exp(e1 - emax) + np.exp(e2 - emax))

    def _lnR2(self, j, h):

        e1 = -2*j
        e2 = 2*j - 4*h
        emax = np.amax([e1, e2])

        return emax + np.log(np.exp(e1 - emax) + np.exp(e2 - emax))

    def _lnR3(self, j, h):

        e1 = 2*h
        e2 = -2*h
        emax = np.amax([e1, e2])

        return emax + np.log(np.exp(e1 - emax) + np.exp(e2 - emax))

    def J(self, interaction, field):
        j = self.m * interaction
        h = self.m * field
        return (1/4) * (self._lnR1(j, h) + self._lnR2(j, h) - 2*self._lnR3(j, h))

    def H(self, interaction, field):
        j = self.m * interaction
        h = self.m * field
        return (1/4) * (self._lnR1(j, h) - self._lnR2(j, h))

    def G(self, interaction, field):
        j = self.m * interaction
        h = self.m * field
        return (1/4) * (self._lnR1(j, h) + self._lnR2(j, h) + 2*self._lnR3(j, h))

    def _recursion_matrix(self, j, h):

        X = self._pderivative(self.G, 0, [j, h])
        Y = self._pderivative(self.J, 0, [j, h])
        Z = self._pderivative(self.H, 0, [j, h])

        A = self._pderivative(self.G, 1, [j, h])
        B = self._pderivative(self.J, 1, [j, h])
        C = self._pderivative(self.H, 1, [j, h])

        return np.array([[self.eigen, X,   A],
                         [0,          Y,   B],
                         [0,          Z,   C]])

    
    def flow(self, interaction, field, n, output=False):

        j, h = interaction, field

        if output:
            output = [[j, h]]
        else:
            print("k   J                H")
            print("______________________")
            print(f"{0}   {j}         {h}")

        j, h = self.J_1(j,h), self.H_1(j,h)
                    
        for i in range(1, n + 1):

            j, h = self.J(j,h), self.H(j,h)

            if output:
                output.append([j, h])
            else:
                print(f"{i}   {j}         {h}")

        if output:
            return np.array(output)

    def critical_point(self, decimal_precision=5):

        p = 0.01
        jc = 0.02

        for n in range(decimal_precision + 1):
            jc = jc - p
            p = 1 / 10**(n + 1)

            while True:
                jx = jc
                jx = self.J(jx, 0)
                if jx < jc:
                    jc = jc + p
                else:
                    break

        return round(jc, decimal_precision)

    def densities(self, field=0):
    # Densities under fixed H

            temp_list = []
            M_results = []

            j_initial = 0.02

            while j_initial < 10:

                h = field
                j = j_initial
                temp_list.append(1 / j)

                U = np.identity(3)
        
                if h == 0:
                    if j < self.Jc:
                        Mn = [1, 0, 0]
                        
                        U = self.neigen * np.dot(self._recursion_matrix_1(j, h), U)
                        j, h = self.J_1(j, h), self.H_1(j, h)
                        while j > 0:
                            U = self.neigen * np.dot(self._recursion_matrix(j, h), U)
                            j, h = self.J(j, h), self.H(j, h)
                        M = np.dot(Mn, U)
                    
                    elif j > self.Jc:
                        Mn = [1, 1, 2]
                        
                        U = self.neigen * np.dot(self._recursion_matrix_1(j, h), U)
                        j, h = self.J_1(j, h), self.H_1(j, h)
                        while j < 1000:
                            U = self.neigen * np.dot(self._recursion_matrix(j, h), U)
                            j, h = self.J(j, h), self.H(j, h)
                        M = np.dot(Mn, U)
                
                else:
                    if h > 0:
                        Mn = [1, 1, 2]
                        
                    elif h < 0:
                        Mn = [1, 1, -2]
                    
                    U = self.neigen * np.dot(self._recursion_matrix_1(j, h), U)
                    j, h = self.J_1(j, h), self.H_1(j, h)
                    
                    js, hs = j, h
                    count = 1
                    while j > 0:
                        j, h = self.J(j, h), self.H(j, h)
                        count = count + 1
                    j, h = js, hs
                    
                    for i in range(count):
                        U = self.neigen * np.dot(self._recursion_matrix(j, h), U)
                        j, h = J(j, h), H(j, h)
                    M = np.dot(Mn, U)
                
                M_results.append(M)

                #j_initial = j_initial + 0.01
                
                if j_initial < self.Jc - 0.01:
                    j_initial = j_initial + 0.001
                
                elif j_initial < self.Jc + 0.01:
                    j_initial = j_initial + 0.0001
                
                elif j_initial < self.Jc + 0.3:
                    j_initial = j_initial + 0.01
                
                elif j_initial < self.Jc + 0.5:
                    j_initial = j_initial + 0.01
                
                else:
                    j_initial = j_initial + 1

            return np.array(temp_list), np.array(M_results)