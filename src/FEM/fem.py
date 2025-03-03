import numpy as np
import utils as MSA

class Node:
    def __init__(self, coords: np.array, BCs: np.ndarray(6) == None, loads: np.zeros(6), id):
        
        self.coords = coords
        self.BCs = BCs
        self.loads = loads
        self.id = id

    def nodal_info(self):
        return self.coords, self.BCs, self.loads, self.id
    
class Element:
    def __init__(self, E, nu, A, Iy, Iz, J, coords1, coords2, id, v_temp):
        self.E = E
        self.nu = nu
        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.J = J
        self.coords1 = coords1
        self.coords2 = coords2
        self.id = id
        self.v_temp = v_temp

    def L(self):
        return np.sqrt((self.coords2[0] - self.coords1[0]) ** 2 + (self.coords2[1] - self.coords1[1]) ** 2 + (self.coords2[2] - self.coords1[2]) ** 2)

    
    def local_K(self):
        return MSA.local_elastic_stiffness_matrix_3D_beam(self.E, self.nu, self.A, self.L(), self.Iy, self.Iz, self.J)

    def gamma(self):
        
        Gamma = MSA.rotation_matrix_3D(self.coords1[0], self.coords1[1], self.coords1[2], self.coords2[0], self.coords2[1], self.coords2[2],self.v_temp)
        return MSA.transformation_matrix_3D(Gamma)

    def global_K(self):
        return self.gamma().T @ self.local_K() @ self.gamma()

    def el_info(self):
        return self.global_K(), self.id
    
class Fem:
    def __init__(self, num_nodes, dof_per_node, K_el, bc, load, id):
        self.num_nodes = num_nodes
        self.dof_per_node = dof_per_node
        self.K_el = K_el
        self.bc = bc
        self.load = load
        self.id = id

    def connectivity(self):
        con  = []
        for i in range(len(self.id)):
            con.append(self.id[i])
        return np.array(con)

    def Big_K(self):
        connectivity = self.connectivity()
        K = np.zeros([self.num_nodes * self.dof_per_node, self.num_nodes * self.dof_per_node])
        for elem in range(connectivity.shape[0]):  # Loop over elements
            nodes = connectivity[elem]  # Get nodes for this element
            Kel = self.K_el[elem]
    
    # Compute global DOF indices for this element
            global_dof_indices = np.concatenate([
                np.arange(6 * nodes[0], 6 * nodes[0] + 6),  # DOFs of first node
                np.arange(6 * nodes[1], 6 * nodes[1] + 6)   # DOFs of second node
            ])
    
    # Assemble local stiffness matrix into the global stiffness matrix
            for p in range(12):  # Local row
                for q in range(12):  # Local column
                    global_p = int(global_dof_indices[p])  # Map local to global row
                    global_q = int(global_dof_indices[q]) # Map local to global column
            
                    K[global_p, global_q] += Kel[p, q]  # Assemble
        return K

    def BC_vec(self):
        vec = np.zeros(self.num_nodes * self.dof_per_node)
        for i in range(self.num_nodes):
            vec[6*i : 6*i + 6] = self.bc[i]
        return np.array([bool(x) for x in vec])

    def force_vec(self):
        vec = np.zeros(self.num_nodes * self.dof_per_node)
        for i in range(self.num_nodes):
            vec[6*i : 6*i + 6] = self.load[i]
        return vec
    def fem_info(self):
        return self.Big_K(), self.BC_vec(), self.force_vec()
    
class solver:
    def __init__(self, num_nodes, dof_per_node, K, BC, F):
        self.num_nodes = num_nodes
        self.dof_per_node = dof_per_node
        self.K = K
        self.BC = BC
        self.F = F

    def solve(self):
        u_total = np.zeros(self.num_nodes * self.dof_per_node)
        id = np.where(self.BC==True)[0]
        all_indices = np.arange(self.num_nodes * self.dof_per_node)
        missing_id = np.setdiff1d(all_indices, id)
        #deleting

        modified_k = np.delete(self.K, id, axis=0)
        modified_k = np.delete(modified_k, id, axis=1)
        modified_f = np.delete(self.F, id)

        u = np.linalg.solve(modified_k, modified_f)
        u_total[missing_id] = u

        forces = self.K @ u_total
        return u_total, forces