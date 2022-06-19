import time

import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
from bc import BoundaryConditions
from d2q9 import D2Q9
from domain_decompose import DomainDecompose
from sim_params import SimParams
from communicate import Communicate
from utils import error
nx = 64
ny = 64
num_halos = 1

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dc = DomainDecompose(nx, ny, comm)
dc.decompose_x(num_halos)

print(rank, dc.tot_num_subdomain_pts_x, dc.tot_num_subdomain_pts_y)
print(rank, dc.global_str_pt_x, dc.global_end_pt_x, dc.global_str_pt_y, dc.global_end_pt_y)

tot_sum = np.array([0])
comm.Reduce([np.array(dc.tot_num_subdomain_pts_x), MPI.INT], [tot_sum, MPI.INT], op=MPI.SUM, root=0)
if rank == 0:
    print(tot_sum)

ndim = 2
ref_len = 1.0
model = D2Q9()
sim_params = SimParams(mach=0.1, reynolds_num=100, ref_len=ref_len, res=ref_len / nx)
q = model.q

rho = 1.0
u = 0
v = 0

f = np.zeros((q, dc.tot_num_subdomain_pts_x, dc.tot_num_subdomain_pts_y))

moments = np.zeros((3, dc.tot_num_subdomain_pts_x, dc.tot_num_subdomain_pts_y))
moments[0] = rho
moments[1] = u
moments[2] = v
u0 = np.zeros((dc.tot_num_subdomain_pts_x, dc.tot_num_subdomain_pts_y))
v0 = np.zeros((dc.tot_num_subdomain_pts_x, dc.tot_num_subdomain_pts_y))

f = model.feq(moments)
moments_old = np.zeros((3, dc.tot_num_subdomain_pts_x, dc.tot_num_subdomain_pts_y))
iter = 0
err = 1.0
start = MPI.Wtime()
# while iter < int(sim_time):
# while err > 1e-6:
while iter < int(20000):
    iter += 1
    # print(iter)

    model.collision(f, sim_params.beta)
    bc = BoundaryConditions(f, comm)
    com = Communicate(f, comm)
    com.communicate_x()
    model.advection(f)

    if rank == 0:
        bc.bounce_back_left()
    if rank == size - 1:
        bc.bounce_back_right()

    bc.bounce_back_bottom()
    bc.moving_wall_bc_top(sim_params.ref_vel)

    if iter % 1000 == 0:
        err = error(f, mom_old=moments_old, comm=comm)
        if rank == 0:
            print("iter: ", iter, "total error: ", err)

        moments_old = model.moments(f)


if rank == 0:
    print("total time taken: ", MPI.Wtime() - start)
# comm.Disconnect()

moments = model.moments(f)
# to plot contours of normal velocity
x = np.linspace(0, 1, dc.tot_num_subdomain_pts_x)
y = np.linspace(0, 1, dc.tot_num_subdomain_pts_y)
X, Y = np.meshgrid(x, y)
VelNorm = np.sqrt(moments[1] * moments[1]
                  + moments[2] * moments[2])

plt.figure()
# cp = plt.contourf(X, Y, np.transpose(VelNorm), 25, cmap=plt.cm.viridis)
cp = plt.contour(X, Y, np.transpose(VelNorm), 15)
cbar = plt.colorbar(cp)
cbar.set_label('Velocity', fontsize=12)
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('velocity_cont_'+str(rank)+'.png', bbox_inches='tight')
plt.title('Velocity', fontsize=12)
MPI.Finalize()


