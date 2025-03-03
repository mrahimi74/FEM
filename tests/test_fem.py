import pytest
import numpy as np
from FEM import fem


def test_node():
    coords = np.array([1.0, 2.0, 3.0])
    BCs = np.array([1, 0, 0, 1, 0, 1])
    loads = np.array([0, 0, -10, 0, 0, 0])
    node = fem.Node(coords, BCs, loads, id=1)
    assert fem.node.nodal_info() == (coords, BCs, loads, 1)

def test_element():
    E, nu, A, Iy, Iz, J = 210e9, 0.3, 0.02, 1e-6, 1e-6, 5e-6
    coords1, coords2 = np.array([0, 0, 0]), np.array([1, 1, 1])
    v_temp = np.array([0, 1, 0])
    elem = fem.Element(E, nu, A, Iy, Iz, J, coords1, coords2, id=1, v_temp=v_temp)
    
    assert np.isclose(elem.L(), np.sqrt(3))
    assert isinstance(elem.local_K(), np.ndarray)
    assert isinstance(elem.global_K(), np.ndarray)
    assert elem.el_info()[1] == 1

def test_fem():
    num_nodes, dof_per_node = 2, 6
    K_el = [np.eye(12)]
    bc = [np.array([1, 1, 0, 0, 0, 0]), np.array([0, 0, 1, 1, 0, 0])]
    load = [np.zeros(6), np.array([0, 0, -10, 0, 0, 0])]
    id = [[0, 1]]
    fem_model = fem.Fem(num_nodes, dof_per_node, K_el, bc, load, id)
    
    assert np.array_equal(fem_model.connectivity(), np.array(id))
    assert isinstance(fem_model.Big_K(), np.ndarray)
    assert isinstance(fem_model.BC_vec(), np.ndarray)
    assert isinstance(fem_model.force_vec(), np.ndarray)

def test_solver():
    num_nodes, dof_per_node = 2, 6
    K = np.eye(12) * 1000
    BC = np.array([1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0], dtype=bool)
    F = np.array([0, 0, -10, 0, 0, 0, 0, 0, -10, 0, 0, 0])
    
    fem_solver = fem.solver(num_nodes, dof_per_node, K, BC, F)
    u_total, forces = fem_solver.solve()
    
    assert isinstance(u_total, np.ndarray)
    assert isinstance(forces, np.ndarray)
    assert u_total.shape == (num_nodes * dof_per_node,)
    assert forces.shape == (num_nodes * dof_per_node,)