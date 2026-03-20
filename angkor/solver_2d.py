# solver_2d.py
# ANGKOR - 2D Two-Group Neutron Diffusion Solver
# Finite difference method with power iteration

# ============================================================
# IMPORTS
# ============================================================
import numpy as np
from scipy import sparse
from scipy.sparse import linalg

# ---- PART 1 - Solve2D Class 
class Solver2D:
    """
    2D two-group neutron diffusion solver. 
    Using finite difference method and power iterations.
    """
    
    def __init__(self, engine, materials, settings, boundary=None):
        self.engine     = engine
        self.materials  = materials 
        self.settings   = settings 
        self.bc = boundary or {
                "left": "vacuum", "right": "vacuum",
                "top": "vacuum",  "bottom": "vacuum"
            }
    
        # Mesh dimension - get from engine 
        self.nx     = len(self.engine.x_centers)
        self.ny     = len(self.engine.y_centers)
        self.N      = self.nx * self.ny     # Total cell per group 
        
        self.dx     = self.engine.domain_x/self.nx 
        self.dy     = self.engine.domain_y/self.ny 
        
        # Result - filled after solve()
        self.k_eff  = None 
        self.flux1  = None # group 1 flux (2D array)
        self.flux2  = None # group 2 flux (2D array)
        
        print(f" Solve2D initialized")
        print(f" Mesh   : {self.nx} * {self.ny} = {self.N} cells")
        print(f" dx     : {self.dy:.4f} cm ")
        print(f" dy     : {self.dy:.4f} cm ")
        
    def _get_xs(self, i, j):
        """Return cross section at cell (i,j) for both groups 

        Args:
            i (_type_): _description_
            j (_type_): _description_
        """
        
        #Step 1: find material name at this cell 
        mat_name = self.engine.material_map[j][i]
        
        #Step 2: get material properties 
        mat = self.materials[mat_name]
        return {
            "D1"            : mat["D1"],
            "D2"            : mat["D2"],
            "sigma_a1"      : mat["sigma_a1"],
            "sigma_a2"      : mat["sigma_a2"],
            "sigma_s12"     : mat["sigma_s12"],
            "nu_sigma_f1"   : mat["nu_sigma_f1"],
            "nu_sigma_f2"   : mat["nu_sigma_f2"]
        }
    
    def _boundary_coeff(self, D, h, side):
            """
            Return boundary leakage coefficient.
            Simple vacuum:      D / h²
            Extrapolated vacuum: D / (h × (h/2 + 2.1312×D))
            Reflective:         0 (no leakage)
            """
            bc_type = self.bc.get(side, "vacuum")
            # print(f"    DEBUG: side ={side}, bc_type={bc_type}, D={D}")
            if bc_type == "reflective":
                return 0.0                               
            elif bc_type == "vacuum":
                # Simple zero-flux BC:
                # return D / h*
                # Extrapolated BC (more accurate):
                d_ex = 2.1312 * D
                return D / (h * (h/2 + d_ex))           
            else:
                return D / h**2                          

    def _build_matrices(self):
        """
        Build spare matrices A (destruction) and F (fission).
        Uses 5-points finite difference stencil.
        
        Matrix size: (2Nx2N) where N=nx*ny 
        """
        N   = self.N 
        nx  = self.nx 
        ny  = self.ny 
        dx  = self.dx 
        dy  = self.dy 
        
        # --- Sparse Matrix Builders ---
        # --- We collect (row, col, value) triplets first 
        # --- Then build sparse matrix at the end 
        rows_A = []         # row indices for matrix A 
        cols_A = []         # column indices for matrix A 
        vals_A = []         # values of matrix A 
        
        rows_F = []         # row indices for matrix F 
        cols_F = []         # column indices for matrix F  
        vals_F = []         # values of matrix F 
        
        # Helper: Add one entry to matrix A 
        def add_A(row, col, val):
            rows_A.append(row)
            cols_A.append(col)
            vals_A.append(val)
            
        # Helper: Add one entry to matrix F
        def add_F(row, col, val):
            rows_F.append(row)
            cols_F.append(col)
            vals_F.append(val)
            
        # --------------------------------------
        # Main Loop - loop over every cell (i,j)
        # --------------------------------------
            
        for j in range(ny):
            for i in range(nx):
                
                # Cell index 
                n = j*nx+i 
                
                # get cross section at this cell 
                xs = self._get_xs(i,j)
                
                # Read cross section 
                D1          = xs["D1"]
                D2          = xs["D2"]
                siga1       = xs["sigma_a1"]
                siga2       = xs["sigma_a2"]
                sigs12      = xs["sigma_s12"]
                nusigf1     = xs["nu_sigma_f1"]
                nusigf2     = xs["nu_sigma_f2"]
                
                # Diffusion coefficient 
                Ex1 = D1/dx**2
                Ey1 = D1/dy**2
                Ex2 = D2/dx**2
                Ey2 = D2/dy**2
                
                # Center coeffiicent 
                C1 = siga1 + sigs12 
                C2 = siga2
                
                # ---------------------------------------
                # Group 1 row: index = n 
                # Group 2 row: index = n + N 
                # ---------------------------------------
                
                # --- Scattering: group 1 to group 2 ----
                # - Scattering Removes from group 1 
                # - Scattering Adds to group 2 equation (negative off-diagonal)
                # ---------------------------------------
                add_A(n+N, n, -sigs12)
                
                # --- Fission source terms --------------
                add_F(n, n, nusigf1)
                add_F(n, n+N, nusigf2) 
                
                # --- Neighbor terms (5 points stencil) -
                # --- West boundary ---
                if i > 0:
                    n_left = j*nx + (i-1)
                    add_A(n,   n_left,   -Ex1)
                    add_A(n+N, n_left+N, -Ex2)
                    C1 += Ex1
                    C2 += Ex2
                else:
                    C1 += self._boundary_coeff(D1, dx, "left")   
                    C2 += self._boundary_coeff(D2, dx, "left")   

                # --- East boundary ---
                if i < nx-1:
                    n_right = j*nx + (i+1)
                    add_A(n,   n_right,   -Ex1)
                    add_A(n+N, n_right+N, -Ex2)
                    C1 += Ex1
                    C2 += Ex2
                else: 
                    C1 += self._boundary_coeff(D1, dx, "right")  
                    C2 += self._boundary_coeff(D2, dx, "right")  

                # --- South boundary ---
                if j > 0:
                    n_down = (j-1)*nx + i
                    add_A(n,   n_down,   -Ey1)
                    add_A(n+N, n_down+N, -Ey2)
                    C1 += Ey1
                    C2 += Ey2
                else:
                    C1 += self._boundary_coeff(D1, dy, "bottom") # ← NEW
                    C2 += self._boundary_coeff(D2, dy, "bottom") # ← NEW

                # --- North boundary ---
                if j < ny-1:
                    n_up = (j+1)*nx + i
                    add_A(n,   n_up,   -Ey1)
                    add_A(n+N, n_up+N, -Ey2)
                    C1 += Ey1
                    C2 += Ey2
                else:
                    C1 += self._boundary_coeff(D1, dy, "top")    # ← NEW
                    C2 += self._boundary_coeff(D2, dy, "top")    # ← NEW
                
                # --- Diagonal (center) terms -----------
                add_A(n, n, C1)
                add_A(n+N, n+N, C2)
    
        # -----------------------------------------
        # --- Build sparse matrix from triplets ---
        # -----------------------------------------
        size = 2*N 
        self.A = sparse.csr_matrix(
            (vals_A, (rows_A, cols_A)),
            shape=(size,size)
        )
        
        self.F = sparse.csr_matrix(
            (vals_F,(rows_F, cols_F)),
            shape=(size,size)
        )
        
        # print(f"  Diagonal sample (boundary cell): {self.A[0,0]:.6f}")
        # print(f"  Diagonal sample (interior cell): {self.A[1000,1000]:.6f}")
        
        # DEBUG — remove after confirming:
        # print(f"  A[0,0] = {self.A[0,0]:.6f}")      # boundary cell
        # print(f"  A[500,500] = {self.A[500,500]:.6f}")  # interior cell
        
        print(f" Matrix A: {self.A.shape}," f"{self.A.nnz} non-zeros")
        print(f" Matrix F: {self.F.shape}," f"{self.F.nnz} non-zeros")
        
        
    def solve(self):
        """
        Run power iteration to find k-eff and flux.
        Algorithm:
        1. Guess k=1.0, flux=ones
        2. Solve A·phi_new = (1/k)·F·phi_old
        3. Update k
        4. Check convergence
        5. Repeat until converged
        """
        
        N        = self.N
        max_iter = self.settings.max_iterations
        tol      = self.settings.convergence
        
        print(f"\n  Starting power iteration...")
        print(f"  Max iterations : {max_iter}")
        print(f"  Tolerance      : {tol}")
        
        # Step 1: Initial guess
        # flux vector size = 2N (group1 + group2)
        phi = np.ones(2 * N)     # flat flux guess
        k   = 1.0                # initial k-eff guess
        
        # Step 2: Build matrices
        print(f"\n  Building matrices...")
        self._build_matrices()
        
        # Step 3: Power iteration loop
        for iteration in range(max_iter):
            # --- Compute fission source: F·phi ---
            fission_source = self.F.dot(phi)
            # --- Solve linear system: A·phi_new = (1/k)·fission_source ---
            rhs     = fission_source / k
            phi_new = linalg.spsolve(self.A, rhs)
            # --- Compute new fission source ---
            fission_new = self.F.dot(phi_new)
            # --- Update k-eff ---
            k_new = k * (np.sum(fission_new) / np.sum(fission_source))
            # --- Check convergence ---
            k_error   = abs(k_new - k)
            phi_error = np.max(np.abs(phi_new - phi) / (np.abs(phi_new) + 1e-10))
            # --- Print progress every 10 iterations ---
            if (iteration + 1) % 10 == 0 or iteration == 0:
                print(f"  Iter {iteration+1:4d}: "
                      f"k={k_new:.6f}  "
                      f"dk={k_error:.2e}  "
                      f"dphi={phi_error:.2e}")
            # --- Update for next iteration ---
            phi = phi_new / np.max(phi_new)   # normalize
            k   = k_new
            # --- Converged? ---
            if k_error < tol and phi_error < tol:
                print(f"\n  CONVERGED at iteration {iteration+1}!")
                break
        # Step 4: Store results
        self.k_eff = k
        # Reshape flux from 1D vector → 2D arrays
        self.flux1 = phi[:N].reshape(self.ny, self.nx)
        self.flux2 = phi[N:].reshape(self.ny, self.nx)
        print(f"\n  {'='*40}")
        print(f"  k-eff = {self.k_eff:.6f}")
        print(f"  {'='*40}")
        return self.k_eff, self.flux1, self.flux2    
                 
    
# TEST 
if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from input_reader import InputReader

    print("=" * 50)
    print("  ANGKOR — 2D Diffusion Solver Test")
    print("=" * 50)

    # Load input
    reader = InputReader("input/pwr_2d.yaml")
    reader.read()

    # Create and run solver
    solver = Solver2D(reader.engine,
                      reader.materials,
                      reader.solver)

    k_eff, flux1, flux2 = solver.solve()

    # Plot flux
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    im1 = axes[0].imshow(flux1, origin="lower",
                          extent=[0, reader.engine.domain_x,
                                  0, reader.engine.domain_y],
                          cmap="hot", aspect="equal")
    axes[0].set_title(f"Group 1 Flux (Fast)\nk-eff = {k_eff:.6f}")
    axes[0].set_xlabel("x (cm)")
    axes[0].set_ylabel("y (cm)")
    plt.colorbar(im1, ax=axes[0])

    im2 = axes[1].imshow(flux2, origin="lower",
                          extent=[0, reader.engine.domain_x,
                                  0, reader.engine.domain_y],
                          cmap="hot", aspect="equal")
    axes[1].set_title(f"Group 2 Flux (Thermal)\nk-eff = {k_eff:.6f}")
    axes[1].set_xlabel("x (cm)")
    axes[1].set_ylabel("y (cm)")
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    plt.savefig("flux_map.png", dpi=150)
    plt.show()
    print("\n  Flux map saved to flux_map.png")
