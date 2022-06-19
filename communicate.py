from mpi4py import MPI
import numpy as np


class Communicate:
    def __init__(self, f: np.ndarray, comm: MPI.COMM_WORLD):
        self.comm = comm
        self.f = f
        self.size = comm.Get_size()
        self.rank = comm.Get_rank()

    def communicate_x(self):
        if self.rank % 2 == 0:
            if self.rank != self.size - 1:
                self.comm.send(self.f[:, -2, :], dest=self.rank + 1, tag=0)
        else:
            self.f[:, 0, :] = self.comm.recv(source=self.rank - 1, tag=0)

        if self.rank % 2 == 1:
            if self.rank != self.size - 1:
                self.comm.send(self.f[:, -2, :], dest=self.rank + 1, tag=1)
        else:
            if self.rank != 0:
                self.f[:, 0, :] = self.comm.recv(source=self.rank - 1, tag=1)

        if self.rank % 2 == 0:
            if self.rank != 0:
                self.comm.send(self.f[:, 1, :], dest=self.rank - 1, tag=2)
        else:
            if self.rank != self.size - 1:
                self.f[:, -1, :] = self.comm.recv(source=self.rank + 1, tag=2)

        if self.rank % 2 == 1:
            self.comm.send(self.f[:, 1, :], dest=self.rank - 1, tag=3)
        else:
            if self.rank != self.size - 1:
                self.f[:, -1, :] = self.comm.recv(source=self.rank + 1, tag=3)
