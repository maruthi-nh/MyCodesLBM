#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from d2q9 import D2Q9


class BoundaryConditions(D2Q9):

    def __init__(self, f):
        self.f = f
        self.f[:, 0, :] = self.f[:, 1, :]
        self.f[:, -1, :] = self.f[:, -2, :]
        self.f[:, :, 0] = self.f[:, :, 1]
        self.f[:, :, -1] = self.f[:, :, -2]
        D2Q9.__init__(self)
        
    def periodic_in_x(self):
        """
        To apply periodic boundary condition in x
        """
        self.f[1, 0, 1:-1] = self.f[1, -1, 1:-1]
        self.f[8, 0, 1:-1] = self.f[8, -1, 1:-1]
        self.f[5, 0, 1:-1] = self.f[5, -1, 1:-1]

        self.f[6, -1, 1:-1] = self.f[6, 0, 1:-1]
        self.f[3, -1, 1:-1] = self.f[3, 0, 1:-1]
        self.f[7, -1, 1:-1] = self.f[7, 0, 1:-1]

    def bounce_back_right(self):
        """
        Apply bounce back on right wall
        """
        self.f[3, -2, 2:-2] = self.f[1, -1, 2:-2]
        self.f[6, -2, 2:-2] = self.f[8, -1, 2:-2]
        self.f[7, -2, 2:-2] = self.f[5, -1, 2:-2]

    def bounce_back_left(self):
        """
        Apply bounce back on right wall
        """
        self.f[8, 1, 2:-2] = self.f[6, 0, 2:-2]
        self.f[1, 1, 2:-2] = self.f[3, 0, 2:-2]
        self.f[5, 1, 2:-2] = self.f[7, 0, 2:-2]

    def bounce_back_bottom(self):
        """
        Apply bounce back on right wall
        """
        self.f[2, 2:-2, 1] = self.f[4, 2:-2, 0]
        self.f[6, 2:-2, 1] = self.f[8, 2:-2, 0]
        self.f[5, 2:-2, 1] = self.f[7, 2:-2, 0]

    def moving_wall_bc_top(self, u_wall):
        """
        Apply moving wall bounce back on top wall
        """
        density=self.f.sum(axis=0)
        self.f[4, 1:-1, -2] = self.f[2, 1:-1, -1]
        self.f[7, 1:-1, -2] = self.f[5, 1:-1, -1] + 6.0 * density[1:-1, -1] * self.w[7] * self.c[0][7] * u_wall
        self.f[8, 1:-1, -2] = self.f[6, 1:-1, -1] + 6.0 * density[1:-1, -1] * self.w[8] * self.c[0][8] * u_wall

    def bot_left_corner_correction(self):
        """
        Bottom left corner
        """
        self.f[2, 1, 1] = self.f[4, 1, 0]
        self.f[5, 1, 1] = self.f[7, 1, 0]
        self.f[1, 1, 1] = self.f[3, 1, 0]
        self.f[8, 1, 1] = self.f[6, 1, 0]
        self.f[6, 1, 1] = self.f[8, 1, 0]

    def bot_right_corner_correction(self):
        """
        Bottom right corner
        """
        self.f[2, -2, 1] = self.f[4, -2, 0]
        self.f[6, -2, 1] = self.f[8, -2, 0]
        self.f[3, -2, 1] = self.f[1, -2, 0]
        self.f[5, -2, 1] = self.f[7, -2, 0]
        self.f[7, -2, 1] = self.f[5, -2, 0]

    def top_left_corner_correction(self):
        """
        Top left corner
        """
        self.f[5, 1, -2] = self.f[7, 1, -1]
        self.f[7, 1, -2] = self.f[5, 1, -1]
        self.f[1, 1, -2] = self.f[3, 1, -1]
        self.f[8, 1, -2] = self.f[6, 1, -1]
        self.f[4, 1, -2] = self.f[2, 1, -1]

# def MovingWallBcOnGridBounceBack(Uwall, Weights, PartVel, PrimVars, f):
        """
        Moving wall boundary conditions. Implemented from Q. Zou, X.He: Phy
    sc Fluids, 1997
    """

    # f[Ny][i][4] = f_post[Ny][i][2]
    # f[Ny][i][7] = f_post[Ny][i][5] + 6.0 * rho[Ny][i] * w[7] * cx[7] * Uwall;
    # f[Ny][i][8] = f_post[Ny][i][6] + 6.0 * rho[Ny][i] * w[8] * cx[8] * Uwall;

    # f[4, 1:-1, -2] = f[2, 1:-1, -1]
    # f[7, 1:-1, -2] = f[5, 1:-1, -1] + 6.0 * PrimVars[0, 1:-1, -1] * Weights[7] * PartVel[0][7] * Uwall
    # f[8, 1:-1, -2] = f[6, 1:-1, -1] + 6.0 * PrimVars[0, 1:-1, -1] * Weights[8] * PartVel[0][8] * Uwall

    # #first set like a solid wall
    # Grids[4, 2:-2, -2] = Grids[2, 2:-2, -1]
    # Grids[8, 2:-2, -2] = Grids[6, 2:-2, -1]
    # Grids[7, 2:-2, -2] = Grids[5, 2:-2, -1]

    # #now apply moving wall boundary conditions
    # PrimVarTemp = np.array([1, Uwall, 0])

    # feqTemp = np.zeros(9)
    # D2Q9.CalEqbDistrFun(PartVel, Weights, PrimVarTemp, feqTemp)
    # #print("priting in the BCs feq", feqTemp)

    # denominator =1.0/(feqTemp[7]+feqTemp[4]+feqTemp[8])

    # factor = Grids[2, 2:-2, -1] + Grids[6, 2:-2, -1] + Grids[5, 2:-2, -1]

    # Grids[4, 2:-2, -2] = factor*denominator*feqTemp[4]
    # Grids[8, 2:-2, -2] = factor*denominator*feqTemp[8]
    # Grids[7, 2:-2, -2] = factor*denominator*feqTemp[7]
    # #print(Grid)

    # rhoTemp = 0
    # rhoTemp = sum(Grids[:,1,-2])
    # #print("rhoTemp=",rhoTemp)
    # Grids[:, 1, -2] = feqTemp*rhoTemp

    # #top right corner
    # Grids[8, -2, -2] = Grids[6, -2, -1]
    # Grids[6, -2, -2] = Grids[8, -2, -1]
    # Grids[3, -2, -2] = Grids[1, -2, -1]
    # Grids[7, -2, -2] = Grids[5, -2, -1]
    # Grids[4, -2, -2] = Grids[2, -2, -1]

    # rhoTemp = 0
    # rhoTemp = sum(Grids[:,-2,-2])
    # Grids[:,-2, -2] = feqTemp*rhoTemp


# def bc_periodic_x(Grid):
#     """
#     To apply periodic BC in x-dir
#     """

    # left wall
    # Grid[1, -1, 1:-1] = Grid[1, 0, 1:-1]
    # Grid[8, -1, 1:-1] = Grid[8, 0, 1:-1]
    # Grid[5, -1, 1:-1] = Grid[5, 0, 1:-1]

    # left wall
    # Grid[6, 0, 1:-1] = Grid[6, -1, 1:-1]
    # Grid[3, 0, 1:-1] = Grid[3, -1, 1:-1]
    # Grid[7, 0, 1:-1] = Grid[7, -1, 1:-1]

    # Grid[1, 0, 1:-1] = Grid[1, -1, 1:-1]
    # Grid[8, 0, 1:-1] = Grid[8, -1, 1:-1]
    # Grid[5, 0, 1:-1] = Grid[5, -1, 1:-1]
    # 
    # Grid[6, -1, 1:-1] = Grid[6, 0, 1:-1]
    # Grid[3, -1, 1:-1] = Grid[3, 0, 1:-1]
    # Grid[7, -1, 1:-1] = Grid[7, 0, 1:-1]



