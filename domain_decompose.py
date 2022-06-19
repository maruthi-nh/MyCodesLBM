from mpi4py import MPI


class DomainDecompose:

    def __init__(self, nx: int, ny: int, comm: MPI.COMM_WORLD):
        self.global_str_pt_x = None
        self.global_end_pt_x = None
        self.global_str_pt_y = None
        self.global_end_pt_y = None
        self.tot_num_subdomain_pts_y = None
        self.tot_num_subdomain_pts_x = None
        self.num_subdomain_pts_x = None
        self.num_subdomain_pts_y = None
        self.nx = nx
        self.ny = ny
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

    def decompose_x(self, num_halos: int):
        """
        To distribute number of points in x among procs

        :param num_halos: Number of halo layers required for data exchange
        :return: number of subdomain pts distributed for each rank, start and end index of each subdomain
        """
        self.num_subdomain_pts_x = self.nx // self.size
        if self.rank < self.nx % self.size:
            self.num_subdomain_pts_x += 1

        self.tot_num_subdomain_pts_x = self.num_subdomain_pts_x+2*num_halos
        self.tot_num_subdomain_pts_y = self.ny + 2*num_halos

        self.global_str_pt_y = 0
        self.global_end_pt_y = self.nx

        if self.rank < self.nx % self.size:
            self.global_str_pt_x = self.rank * self.num_subdomain_pts_x
            self.global_end_pt_x = (self.rank + 1) * self.num_subdomain_pts_x
        else:
            self.global_str_pt_x = (self.nx % self.size) + self.rank * self.num_subdomain_pts_x
            self.global_end_pt_x = (self.nx % self.size) + (self.rank + 1) * self.num_subdomain_pts_x

    def decompose_y(self, num_halos: int):
        """
        To distribute number of points in y among procs

        :return: number of subdomain pts distributed for each rank, start and end index of each subdomain
        """
        self.num_subdomain_pts_y = self.ny // self.size
        if self.rank < self.ny % self.size:
            self.num_subdomain_pts_y += 1

        self.tot_num_subdomain_pts_y = self.num_subdomain_pts_y + 2*num_halos
        self.tot_num_subdomain_pts_x = self.nx + 2*num_halos

        self.global_str_pt_x = 0
        self.global_end_pt_x = self.nx

        if self.rank < self.ny % self.size:
            self.global_str_pt_y = self.rank * self.num_subdomain_pts_y
            self.global_end_pt_y = (self.rank + 1) * self.num_subdomain_pts_y
        else:
            self.global_str_pt_y = (self.ny % self.size) + self.rank * self.num_subdomain_pts_y
            self.global_end_pt_y = (self.ny % self.size) + (self.rank + 1) * self.num_subdomain_pts_y


# dc = DomainDecompose(65, 69, MPI.COMM_WORLD)
# dc.decompose_y(1)
# print(dc.rank, dc.num_subdomain_pts_x, dc.tot_num_subdomain_pts_x, dc.global_str_pt, dc.global_end_pt)
# print(dc.rank, dc.num_subdomain_pts_y, dc.tot_num_subdomain_pts_y, dc.global_str_pt, dc.global_end_pt)
