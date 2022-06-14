#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some important formulas:
    1. formula for kinematic viscosity nu = ((2*tau - 1)*(dx^2))/(6*dt)
    2. lattice speed: lambda = dx/dt
    3. equilibrium distribution function is:
        feq_i = w_i*{3*[(e_i.u)/Lambda] + (9/2)*[(e_i.u)^2/Lambda^2] - (3/2)*[(e_i.u)/Lambda^2]}
"""

import numpy as np
import matplotlib.pyplot as plt

from d2q9 import D2Q9, advection
from bc import BoundaryConditions
from sim_params import SimParams

# import data_output_functions as dop


nx = 100
ny = 100
ghost_x = 2
ghost_y = 2
tot_nx = nx + ghost_x
tot_ny = ny + ghost_y
ndim = 2
model = D2Q9()
sim_params = SimParams(mach=0.1, reynolds_num=100, ref_len=1.0, res=0.02)
q = model.q

sim_time = 20000  # (40*Nx)/UtopWall

rho_ini = 1.0
p = rho_ini * model.T0
u_ini = 0.0
v_ini = 0.0

# initialization of arrays required
f = np.zeros((q, tot_nx, tot_ny))
moments = np.zeros((3, tot_nx, tot_ny))
moments[0] = rho_ini
moments[1] = u_ini
moments[2] = v_ini

u0 = np.zeros((tot_nx, tot_ny))
v0 = np.zeros((tot_nx, tot_ny))

######First iteration#####

f = model.feq(moments)
bc = BoundaryConditions(f)
iter = 0
err = 1.0

# while iter < int(sim_time):
while(err>1e-6):

    # while(iter<int(100)):
    iter = iter + 1
    # print(iter)

    model.collision(f, sim_params.beta)
    BoundaryConditions(f)
    advection(f)
    bc.bounce_back_left()
    bc.bounce_back_bottom()
    bc.bounce_back_right()
    bc.moving_wall_bc_top(sim_params.ref_vel)


    if iter % 1000 == 0:
        moments=model.moments(f)
        e1 = np.sum(np.sqrt((moments[1] - u0) ** 2 + (moments[2] - v0) ** 2))
        e2 = np.sum(np.sqrt(moments[1] ** 2 + moments[2] ** 2))
        err = e1 / e2
        print('error is = ', err)
    u0 = moments[1]
    v0 = moments[2]

# to plot contours of normal velocity
x = np.linspace(1, nx, nx)
y = np.linspace(1, ny, ny)
X, Y = np.meshgrid(x, y)
VelNorm = np.sqrt(moments[1, 1:-1, 1:-1] * moments[1, 1:-1, 1:-1] + moments[2, 1:-1, 1:-1] * moments[2, 1:-1, 1:-1])

cp = plt.contourf(X, Y, np.transpose(VelNorm), 10)
plt.colorbar(cp)
plt.xlabel('x')
plt.ylabel('y')

plt.figure()
plt.plot(moments[1, int(nx / 2), 1:-1])
