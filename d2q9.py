#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


class D2Q9:

    def __init__(self):
        self.c = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],
                           [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
        self.w = np.array([4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                           1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0])
        self.q = 9
        self.T0 = 1 / 3.0
        self.a = np.sqrt(self.T0)

    def feq(self, moments):
        """
        To compute distribution function. Takes particle velocity, weights, and Primitive variables
        (rho, u, v) and gives Feq (equilibrium distribution as output. PrimVars and Feq can be 1D array
        or matrix (nD array)
        """

        f = np.zeros((self.q, moments.shape[1], moments.shape[2]))

        for i in range(self.q):
            udotc = self.c[0][i] * moments[1] + self.c[1][i] * moments[2]
            usq = moments[1] * moments[1] + moments[2] * moments[2]

            f[i] = self.w[i] * moments[0] * (1.0 + 3.0 * udotc - 1.5 * usq + 4.5 * udotc ** 2
                                        + 4.5 * udotc ** 3 - 4.5 * usq * udotc)
        return f

    def moments(self, f):
        """
        Get (rho, u, p) by taking moment of distribution function. Takes particle velocity, and Grids (distribution
        fun) as input and gives PrimVars (rho, u, v) as output
        """
        mom = np.zeros((3, f.shape[1], f.shape[2]))

        mom[0] = f.sum(axis=0)

        for i in range(self.q):
            mom[1] += f[i] * self.c[0, i]
            mom[2] += f[i] * self.c[1, i]

        mom[1] /= mom[0]
        mom[2] /= mom[0]

        return mom

    def collision(self, f, beta):
        """
        computes collision term of the LBM update formula, i.e., f*=f+(alpha*beta)*(feq-f). Takes beta, Particle Vel,
        Weights, Grid, PrimVar and Feq as inputs. These gives updated Grid as output
        """
        # GetMoments(PartVel, Grids, PrimVars)

        feq = self.feq(self.moments(f))
        # vp.feq_iso(PartVel, Weights, PrimVars,  Feq)
        # CalEqbDistrFun_Iso(PartVel, Weights, PrimVars, Feq)

        alpha = 2.0
        for i in range(self.q):
            f[i] += alpha * beta * (feq[i] - f[i])
        

def advection(f):
    """
    Computes Grids after streaming
    """
    len0, lenx, leny = np.shape(f)

    for j in range(0, leny - 1):
        for i in range(0, lenx - 1):
            f[3, i, j] = f[3, i + 1, j]
            f[4, i, j] = f[4, i, j + 1]
            f[7, i, j] = f[7, i + 1, j + 1]
            f[8, i, j] = f[8, i - 1, j + 1]

    for j in range(leny - 2, 0, -1):
        for i in range(lenx - 2, 0, -1):
            f[1, i, j] = f[1, i - 1, j]
            f[2, i, j] = f[2, i, j - 1]
            f[5, i, j] = f[5, i - 1, j - 1]
            f[6, i, j] = f[6, i + 1, j - 1]
    

