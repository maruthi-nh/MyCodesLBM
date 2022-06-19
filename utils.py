from mpi4py import MPI
import numpy as np
from d2q9 import D2Q9


def error(f: np.ndarray, mom_old: np.ndarray, comm: MPI.COMM_WORLD) -> float:
    """
    To compute error between old and new solution
    :param f: distribution
    :param mom_old: old moments
    :param comm: MPI communicator
    :return: error
    """
    model = D2Q9()
    moments = model.moments(f)
    e1 = np.sum(np.sqrt((moments[1] - mom_old[1]) ** 2 + (moments[2] - mom_old[2]) ** 2))
    e2 = np.sum(np.sqrt(moments[1] ** 2 + moments[2] ** 2))
    global_err = np.array([1, 1], 'd')
    comm.Reduce([np.array([e1, e2]), MPI.DOUBLE], [global_err, MPI.DOUBLE], op=MPI.SUM, root=0)

    return global_err[0] / global_err[1]
