#- solver_mg.py
#- ANGKOR - Multi-group 2D neutron diffusion solver 
#- Supports any number of groups G 

import numpy as np 
from scipy import sparse
from scipy.sparse import linalg

class SolverMG:
    """
    G-Group 2D neutron diffusion solver.
    Generalizes Solver2D to any number of energy groups. 
    """
    def __init__(self, engine, materials, settings, n_groups, boundary=None):
        self.engine     = engine 
        self.materials  = materials
        self.settings   = settings 
        self.G          = n_groups 
        self.nx         = len(self.engine.x_centers)
        self.ny         = len(self.engine.y_centers)
        self.N          = self.nx * self.ny 
        self.dx         = self.engine.domain_x /self.nx 
        self.dy         = self.engine.domain_y /self.ny 
        self.buckling   = getattr(settings, "buckling", 0.0)
        self.bc         = boundary or {
            "left": "vacuum",
            "right": "vacuum",
            "top": "vacuum",
            "bottom": "vacuum",
        }
        
        # Results
        self.k_eff      = None 
        self.flux       = None 
        
        print(f"    SolverMG Initialized---")
        print(f"    Groups  : {self.G}")
        print(f"    Mesh    : {self.nx} x {self.ny} = {self.N} cells")
    
    def _get_xs(self, i,j):
        """ 
        Return cross section at cell (i,j) for all groups.
        """
        mat_name    = self.engine.material_map[j][i]
        mat         = self.materials[mat_name]
        return mat 
    
    def _cell_index(self, g, i, j):
        """ 
        Convert (group, i, j) to global matrix index.
        In 2-gorup: group 1 = n, group 2 = n+N 
        in G-group: group g = g*N+n 
                    where n = j*nx+i
        """
        n   = j*self.nx + i
        return g*self.N+n 
    
    def _build_matrices(self):
        """ 
        Build sparse A and F matrices for G-group diffusion.
        """
        size    = self.G * self.N       # total matrix size 
        rows_A, cols_A, vals_A = [], [], []   # loss/leakage
        rows_F, cols_F, vals_F = [], [], []   # fission source
        rows_S, cols_S, vals_S = [], [], []   # scatter-in source
        
        def add_A(r,c,v):
            rows_A.append(r)
            cols_A.append(c)
            vals_A.append(v)
            
        def add_F(r,c,v):
            rows_F.append(r)
            cols_F.append(c)
            vals_F.append(v)

        def add_S(r,c,v):
            rows_S.append(r)
            cols_S.append(c)
            vals_S.append(v)
        
        for j in range(self.ny):
            for i in range(self.nx):
                xs  = self._get_xs(i,j)
                D           = xs["D"]
                sigma_a     = xs["sigma_a"]
                nusigf      = xs["nu_sigma_f"]
                sigma_s     = xs["sigma_s"]     
                chi         = xs["chi"] 
                sigma_a_eff = [sigma_a[g] + D[g]*self.buckling for g in range(self.G)]
                
                for g in range(self.G): 
                    row     = self._cell_index(g,i,j)
                    Dx = D[g] / self.dx**2
                    Dy = D[g] / self.dy**2 
                    
                    # Center coefficeint C=leakage+absoption+All down scattering out
                    sigma_r = sigma_a_eff[g] + sum(sigma_s[g][g2] for g2 in range(self.G) if g2 !=g)
                    C = sigma_r 
                    
                    # x - neighbors
                    if i>0:
                        add_A(row, self._cell_index(g, i-1,j),-Dx)
                        C +=Dx 
                    else:  # left boundary
                        if self.bc.get("left","vacuum") == "vacuum":
                            C += Dx
                    if i < self.nx - 1:
                        add_A(row, self._cell_index(g, i+1, j), -Dx)
                        C += Dx
                    else:  # right boundary
                        if self.bc.get("right","vacuum") == "vacuum":
                            C += Dx

                    # y- neighbors
                    if j > 0:
                        add_A(row, self._cell_index(g, i, j-1), -Dy)
                        C += Dy
                    else:  # bottom boundary
                        if self.bc.get("bottom","vacuum") == "vacuum":
                            C += Dy
                    if j < self.ny - 1:
                        add_A(row, self._cell_index(g, i, j+1), -Dy)
                        C += Dy
                    else:  # top boundary
                        if self.bc.get("top","vacuum") == "vacuum":
                            C += Dy

                    # Diagonal
                    add_A(row, row, C)

                    # Scatter-in source stays on RHS
                    for g_from in range(self.G):
                        if g_from == g:
                            continue 
                        add_S(row, self._cell_index(g_from,i,j), sigma_s[g_from][g])

                    # Fission source 
                    for g2 in range(self.G):
                        add_F(row, self._cell_index(g2, i,j), chi[g]*nusigf[g2])
                        
        self.A = sparse.csr_matrix(
                        (vals_A, (rows_A, cols_A)), shape=(size, size))
        self.F = sparse.csr_matrix(
                        (vals_F, (rows_F, cols_F)), shape=(size, size))
        self.S = sparse.csr_matrix(
                        (vals_S, (rows_S, cols_S)), shape=(size, size))
        print(f"  Matrix A: {self.A.shape}, {self.A.nnz} non-zeros")
        print(f"  Matrix F: {self.F.shape}, {self.F.nnz} non-zeros")
        print(f"  Matrix S: {self.S.shape}, {self.S.nnz} non-zeros")
        print("   Diagonal average leakage term:", np.mean(self.A.diagonal()))
    
    
    def solve(self):
        """
        Run power iteration to find k-eff and flux.
        Same algorithm as Solver2D but for G groups.
        """
        N        = self.N
        G        = self.G
        max_iter = getattr(self.settings, "max_iterations", 1000)
        tol      = self.settings.convergence

        print(f"\n  Starting power iteration (G={G} groups)...")

        # Initial guess — flat flux for all groups
        phi = np.ones(G * N)
        k   = 1.0

        # Build matrices
        print(f"  Building matrices...")
        self._build_matrices()

        # Power iteration loop
        for iteration in range(max_iter):
            fission_source  = self.F.dot(phi)
            scatter_source  = self.S.dot(phi)
            rhs             = scatter_source + fission_source / k
            phi_new         = linalg.spsolve(self.A, rhs) 

            # Update k-eff
            fission_new     = self.F.dot(phi_new) 
            k_new           = k * (np.sum(fission_new) / np.sum(fission_source)) 

            # Check convergence
            k_error         = abs(k_new - k)
            phi_error       = np.max(
                np.abs(phi_new - phi) / (np.abs(phi_new) + 1e-12)
            )

            # Print progress
            if (iteration + 1) % 10 == 0 or iteration == 0:
                print(f"  Iter {iteration+1:4d}: "
                      f"k={k_new:.6f}  "
                      f"dk={k_error:.2e}")

            # Update
            phi = phi_new / np.linalg.norm(phi_new)
            k   = k_new

            # Converged?
            if k_error < tol and phi_error < tol:
                print(f"\n  CONVERGED at iteration {iteration+1}!")
                break

        # Store results
        self.k_eff  = k
        self.flux   = np.zeros((G, self.ny, self.nx))
        for g in range(G):
            self.flux[g] = phi[g*N:(g+1)*N].reshape(self.ny, self.nx)
            #                  ↑   ↑
            # BLANK 7: start = g*N
            # BLANK 8: end   = (g+1)*N

        print(f"\n  {'='*40}")
        print(f"  k-eff = {self.k_eff:.6f}")
        print(f"  {'='*40}")
        return self.k_eff, self.flux

if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from input_reader import InputReader

    print("=" * 50)
    print("  ANGKOR — Multi-Group Solver Test")
    print("=" * 50)

    filepath = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "input", "pwr_2d_4g.yaml")
    reader = InputReader(filepath)
    reader.read()

    solver = SolverMG(
        engine    = reader.engine,
        materials = reader.materials,
        settings  = reader.solver,
        n_groups  = 4
    )

    k_eff, flux = solver.solve()

    # Plot all group fluxes
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    fig.suptitle(f"ANGKOR 4-Group Flux  |  k-eff = {k_eff:.6f}")

    for g in range(4):
        im = axes[g].imshow(
            flux[g],
            origin="lower",
            extent=[0, reader.engine.domain_x,
                    0, reader.engine.domain_y],
            cmap="hot", aspect="equal"
        )
        axes[g].set_title(f"Group {g+1}")
        axes[g].set_xlabel("x (cm)")
        plt.colorbar(im, ax=axes[g])

    plt.tight_layout()
    plt.savefig("flux_4group.png", dpi=150)
    plt.show()
    print("  Saved: flux_4group.png")
