import yaml
from geometry_2d import Region, GeometryEngine 

# ============================================
# SolverSettings class
# ============================================
class SolverSettings:
    """stores solver configuration from input file"""
    def __init__(self, max_iterations, convergence, buckling=0.0):
        self.max_iterations = max_iterations
        self.convergence    = convergence
        self.buckling       = buckling
        
    def describe(self):
        print(f" Max Iterations :   {self.max_iterations}")
        print(f" Convergence    :   {self.convergence}")
        print(f" Buckling (Bz^2):   {self.buckling}")
# ===========================================
# InputReader class
# ===========================================
class InputReader:
    """ 
    Reads an ANGKOR YAML input fiel.
    Builds GeometryEngine and MaterialLibrary read for solver.
    """
    
    def __init__(self, filepath):
        self.filepath  = filepath
        self.data       = None 
        self.engine     = None
        self.materials  = {}
        self.solver     = None 
        self.title      = "Untitled"
        self.boundary   = {
            "left": "vacuum",
            "right": "vacuum",
            "top": "vacuum",
            "bottom": "vacuum",
        }
        
    def read(self):
        """ 
        Our main method to read YAML file and build all objects.
        Call YAML fiel first before accessing engine or materials.
        """
        with open(self.filepath, "r", encoding="utf-8") as f:
            self.data = yaml.safe_load(f)
        
        self.title = self.data.get("title", "Untitled")
        
        self._parse_geometry()
        self._parse_materials()
        self._parse_solver()
        self._parse_boundary()
        
        print(f" Input file read successfully!")
        print(f" TITLE: {self.title}")
        
    def _parse_geometry(self):
        """Build GeometryEngine from YAML geometry+region sections."""
        geo = self.data["geometry"]
        domain_x    = geo["domain_x"]
        domain_y    = geo["domain_y"]
        self.nx     = geo["nx"]
        self.ny     = geo["ny"]
        
        self.engine = GeometryEngine(domain_x, domain_y)
        
        for r in self.data["regions"]:
            region = Region(
                name    = r["name"],
                x_min   = r["x_min"],
                x_max   = r["x_max"],
                y_min   = r["y_min"],
                y_max   = r["y_max"],
                material= r["material"]
            )
            self.engine.add_region(region)
            # Debug print for each region
            # print(f"    Added: {r['name']}  x={r['x_min']}→{r['x_max']}  "
            #       f"y={r['y_min']}→{r['y_max']}  mat={r['material']}")
        self.engine.build_mesh(self.nx, self.ny)
        # Final confirmation
        # print(f"  Final regions in engine: {len(self.engine.regions)}")
    
    def _parse_materials(self):
        """Read material cross sections from YAML input file"""
        mats = self.data["materials"] 
            
        for mat_name, props in mats.items():

        # Detect format by checking which keys exist
            if "groups" in props:
                # --- MULTIGROUP FORMAT ---
                self.materials[mat_name] = {
                    "groups"     : props["groups"],
                    "D"          : props["D"],
                    "sigma_a"    : props["sigma_a"],
                    "nu_sigma_f" : props["nu_sigma_f"],
                    "sigma_s"    : props["sigma_s"],
                    "chi"        : props["chi"]
                }
            else:
                # --- 2-GROUP FORMAT ---
                self.materials[mat_name] = {
                    "D1"          : props["D1"],
                    "D2"          : props["D2"],
                    "sigma_a1"    : props["sigma_a1"],
                    "sigma_a2"    : props["sigma_a2"],
                    "sigma_s12"   : props["sigma_s12"],
                    "nu_sigma_f1" : props["nu_sigma_f1"],
                    "nu_sigma_f2" : props["nu_sigma_f2"],
                }
        self.engine
        
    
    def _parse_solver(self):
        """Read solver settings from YAML."""

        sol = self.data.get("solver", {})  # .get() = safe lookup with default

        self.solver = SolverSettings(
            max_iterations = sol.get("max_iterations", 1000),
            convergence    = sol.get("convergence", 1.0e-6),
            buckling       = sol.get("buckling", 0.0)
        )

    def _parse_boundary(self):
        """Read boundary conditions (optional)."""
        bc = self.data.get("boundary_conditions", {})
        for side in ["left", "right", "top", "bottom"]:
            if side in bc:
                self.boundary[side] = bc[side].lower()

    def summary(self):
        """Print a full summary of the loaded input."""
        print(f"\n{'='*50}")
        print(f"  ANGKOR Input Summary")
        print(f"{'='*50}")
        print(f"  Title    : {self.title}")
        print(f"  File     : {self.filepath}")
        print(f"\n  Geometry:")
        self.engine.summary()
        print(f"\n  Materials: {list(self.materials.keys())}")
        print(f"\n  Solver:")
        self.solver.describe()
        print(f"\n  Boundary conditions: {self.boundary}")

if __name__ == "__main__":

    print("=" * 50)
    print("  TESTING: InputReader")
    print("=" * 50)

    reader = InputReader("input/pwr_2d.yaml")
    reader.read()
    reader.summary()

    # Test material lookup
    print("\n--- Material Check ---")
    print(f"  uo2_fuel D1 : {reader.materials['uo2_fuel']['D1']}")  # 1.40
    print(f"  water D1    : {reader.materials['water']['D1']}")     # 1.13
