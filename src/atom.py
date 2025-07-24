import numpy as np
from math_aux import *

class Atom:
    def __init__(self, n2, l2):
        self.n = n2
        self.l = l2

    def E_Rb(self):
        n = self.n
        l = self.l
        muns_Rb = 3.1311804 + 0.1745312 * (n - 3.1311804) ** -2
        Ens_Rb = -0.5 * (n - muns_Rb) ** -2
        munp_Rb = 2.6482793 + 0.2925324 * (n - 2.6482793) ** -2
        Enp_Rb = -0.5 * (n - munp_Rb) ** -2
        mund_Rb = 1.3472787 - 0.5994376 * (n - 1.3472787) ** -2
        End_Rb = -0.5 * (n - mund_Rb) ** -2
        En = -0.5 * n ** -2
        if l == 0:
            return Ens_Rb
        elif l == 1:
            return Enp_Rb
        elif l == 2:
            return End_Rb
        else:
            return En

    def Vfield(self, li, lj, mi, mj, radial, strength):
        return strength * self.Angular_dc_field(li, lj, mi, mj) * radial

    def Angular_dc_field(self, l, l1, m, m1):
        term1 = 0.0
        if m == m1:
            if l == l1 - 1:
                term1 = np.sqrt((l1 + m1) * (l1 - m1) / (4.0 * l1 * l1 - 1))
            if l == l1 + 1:
                term1 = np.sqrt((l1 + m1 + 1) * (l1 - m1 + 1) / ((2.0 * l1 + 1) * (2.0 * l1 + 3)))
            return term1
        else:
            return 0.0 