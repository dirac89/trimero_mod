import numpy as np
from math_aux import Spherical, DRnl, DOlm, DPhilm, hydrogenicR

class FermiPotentials:
    def __init__(self, s, n, n2, li, lj, mi, mj, r1, theta1, As1, wave1, wave2, Ap1, Dwave1, Dwave2):
        self.s_ = s
        self.n_ = n
        self.n2_ = n2
        self.li_ = li
        self.lj_ = lj
        self.mi_ = mi
        self.mj_ = mj
        self.r1_ = r1
        self.theta1_ = theta1
        self.as1_ = As1
        self.wave1_ = wave1
        self.wave2_ = wave2
        self.ap1_ = Ap1
        self.dWave1_ = Dwave1
        self.dWave2_ = Dwave2

    def Vs(self):
        x1 = np.cos(self.theta1_)
        if self.li_ > 2 and self.lj_ > 2:
            Vs1 = 2 * np.pi * self.as1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * hydrogenicR(self.n_, self.li_, 1., self.r1_) * hydrogenicR(self.n2_, self.lj_, 1., self.r1_)
            return Vs1
        if self.li_ <= 2 and self.lj_ > 2:
            Vs1 = 2 * np.pi * self.as1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * hydrogenicR(self.n2_, self.lj_, 1., self.r1_) * self.wave1_
            return Vs1
        if self.li_ > 2 and self.lj_ <= 2:
            Vs1 = 2 * np.pi * self.as1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * hydrogenicR(self.n_, self.li_, 1., self.r1_) * self.wave2_
            return Vs1
        else:
            Vs1 = 2 * np.pi * self.as1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * self.wave1_ * self.wave2_
            return Vs1

    def Vp(self):
        x1 = np.cos(self.theta1_)
        if self.li_ > 2 and self.lj_ > 2:
            VpA1 = 6 * np.pi * self.ap1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * DRnl(self.n_, self.li_, self.r1_) * DRnl(self.n2_, self.lj_, self.r1_)
            VpB1 = 6 * np.pi * self.ap1_ * DOlm(self.li_, self.mi_, self.theta1_) * DOlm(self.lj_, self.mj_, self.theta1_) \
                * hydrogenicR(self.n2_, self.lj_, 1., self.r1_) * hydrogenicR(self.n_, self.li_, 1., self.r1_) * self.r1_ ** -2
            VpC1 = 6 * np.pi * self.ap1_ * DPhilm(self.li_, self.mi_, self.theta1_) * DPhilm(self.lj_, self.mj_, self.theta1_) \
                * hydrogenicR(self.n_, self.li_, 1., self.r1_) * hydrogenicR(self.n2_, self.lj_, 1., self.r1_) * self.r1_ ** -2
            return VpA1 + VpB1 + VpC1
        if self.li_ <= 2 and self.lj_ > 2:
            VpA1 = 6 * np.pi * self.ap1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * DRnl(self.n2_, self.lj_, self.r1_) * self.dWave1_
            VpB1 = 6 * np.pi * self.ap1_ * DOlm(self.li_, self.mi_, self.theta1_) * DOlm(self.lj_, self.mj_, self.theta1_) \
                * hydrogenicR(self.n2_, self.lj_, 1., self.r1_) * self.wave1_ * self.r1_ ** -2
            VpC1 = 6 * np.pi * self.ap1_ * DPhilm(self.li_, self.mi_, self.theta1_) * DPhilm(self.lj_, self.mj_, self.theta1_) \
                * hydrogenicR(self.n2_, self.lj_, 1., self.r1_) * self.wave1_ * self.r1_ ** -2
            return VpA1 + VpB1 + VpC1
        if self.li_ > 2 and self.lj_ <= 2:
            VpA1 = 6 * np.pi * self.ap1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * DRnl(self.n_, self.li_, self.r1_) * self.dWave2_
            VpB1 = 6 * np.pi * self.ap1_ * DOlm(self.li_, self.mi_, self.theta1_) * DOlm(self.lj_, self.mj_, self.theta1_) \
                * hydrogenicR(self.n_, self.li_, 1., self.r1_) * self.wave2_ * self.r1_ ** -2
            VpC1 = 6 * np.pi * self.ap1_ * DPhilm(self.li_, self.mi_, self.theta1_) * DPhilm(self.lj_, self.mj_, self.theta1_) \
                * hydrogenicR(self.n_, self.li_, 1., self.r1_) * self.wave2_ * self.r1_ ** -2
            return VpA1 + VpB1 + VpC1
        else:
            VpA1 = 6 * np.pi * self.ap1_ * Spherical(self.li_, self.mi_, x1) * Spherical(self.lj_, self.mj_, x1) \
                * self.dWave1_ * self.dWave2_
            VpB1 = 6 * np.pi * self.ap1_ * DOlm(self.li_, self.mi_, self.theta1_) * DOlm(self.lj_, self.mj_, self.theta1_) \
                * self.wave1_ * self.wave2_ * self.r1_ ** -2
            VpC1 = 6 * np.pi * self.ap1_ * DPhilm(self.li_, self.mi_, self.theta1_) * DPhilm(self.lj_, self.mj_, self.theta1_) \
                * self.wave1_ * self.wave2_ * self.r1_ ** -2
            return VpA1 + VpB1 + VpC1

    def Vsp(self):
        return self.Vs() + self.s_ * self.Vp() 