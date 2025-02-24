import pytest
from FEM import fem
import numpy as np


def test_local_elastic_stiffness_matrix_3D_beam():
    E, nu, A, L, Iy, Iz, J = 200e9, 0.3, 0.01, 2.0, 1e-6, 1e-6, 5e-6
    k_e = fem.local_elastic_stiffness_matrix_3D_beam(E, nu, A, L, Iy, Iz, J)
    assert k_e.shape == (12, 12)
    assert np.allclose(k_e, k_e.T), "Matrix should be symmetric"
    assert np.all(np.diag(k_e) >= 0), "Diagonal elements should be non-negative"

def test_check_unit_vector():
    fem.check_unit_vector(np.array([1, 0, 0]))  # Should not raise an error
    with pytest.raises(ValueError):
        fem.check_unit_vector(np.array([2, 0, 0]))

def test_check_parallel():
    vec1 = np.array([1, 0, 0])
    vec2 = np.array([0, 1, 0])
    fem.check_parallel(vec1, vec2)  # Should not raise an error
    with pytest.raises(ValueError):
        fem.check_parallel(vec1, np.array([2, 0, 0]))

def test_rotation_matrix_3D():
    gamma = fem.rotation_matrix_3D(0, 0, 0, 1, 0, 0)
    assert gamma.shape == (3, 3)
    assert np.allclose(np.linalg.det(gamma), 1), "Rotation matrix should have determinant 1"
    assert np.allclose(gamma @ gamma.T, np.eye(3)), "Rotation matrix should be orthogonal"

def test_transformation_matrix_3D():
    gamma = np.eye(3)
    Gamma = fem.transformation_matrix_3D(gamma)
    assert Gamma.shape == (12, 12)
    assert np.allclose(Gamma[0:3, 0:3], gamma)

def test_local_geometric_stiffness_matrix_3D_beam():
    L, A, I_rho, Fx2, Mx2, My1, Mz1, My2, Mz2 = 2.0, 0.01, 5e-6, 100, 10, 5, 5, 5, 5
    k_g = fem.local_geometric_stiffness_matrix_3D_beam(L, A, I_rho, Fx2, Mx2, My1, Mz1, My2, Mz2)
    assert k_g.shape == (12, 12)
    assert np.allclose(k_g, k_g.T), "Matrix should be symmetric"

def test_local_geometric_stiffness_matrix_3D_beam_without_interaction_terms():
    L, A, I_rho, Fx2 = 2.0, 0.01, 5e-6, 100
    k_g = fem.local_geometric_stiffness_matrix_3D_beam_without_interaction_terms(L, A, I_rho, Fx2)
    assert k_g.shape == (12, 12)
    assert np.allclose(k_g, k_g.T), "Matrix should be symmetric"

def test_K_assembly():
    num_nodes = 3
    dof_per_node = 6
    k_element = np.eye(12)  # Simple identity matrix as placeholder
    
    K = fem.K_assembly(num_nodes, dof_per_node, k_element)
    expected_size = num_nodes * dof_per_node
    
    assert K.shape == (expected_size, expected_size)
    assert np.allclose(K, K.T), "Global stiffness matrix should be symmetric"