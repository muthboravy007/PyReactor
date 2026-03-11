# geometry_2d.py
# PyReactor - 2D Rectangular Geometry Engine
# Author: MUTH Boravy
# Description: Defines 2D rectangular regions, builds mesh,
#              validates geometry, and visualizes the core layout

# ============================================================
# IMPORTS
# ============================================================

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.patches as patches 
import matplotlib.colors as mcolors 

# ============================================================
# PART 1: Region Class
# ============================================================
class Region:
    """
    A rectangular region in 2D space with an assigned material.
    Define by its boundary coordinates (x_min, x_max, y_min, y_max)
    """
    def __init__(self, name, x_min, x_max, 
                 y_min, y_max, material):
        self.name       = name 
        self.x_min      = x_min 
        self.x_max      = x_max
        self.y_min      = y_min 
        self.y_max      = y_max 
        self.material   = material
        self._validate()

    def _validate(self):
        """Raise error if region dimension are invalid."""
        if self.x_min >= self.x_max: 
            raise ValueError(
                f"Region '{self.name}': "
                f"x_min ({self.x_min}) must be < x_max ({self.x_max})"
            )
        if self.y_min >= self.y_max: 
            raise ValueError(
                f"Region '{self.name}': "
                f"y_min ({self.y_min}) must be < y_max ({self.y_max})"
            )
    
    
    def contains(self,x,y):
        x_inside = self.x_min <= x <= self.x_max 
        y_inside = self.y_min <= y <= self.y_max 
        return x_inside and y_inside 
    
    def width(self):
        'Return width of region in cm'
        return self.x_max - self.x_min 
    
    def height(self):
        'Return height of region in cm.'
        return self.y_max - self.y_min 
    
    def area(self):
        'Return area of the region in cm^2.'
        return self.width()*self.height() 
    
    def describe(self):
        print(f"Region      : {self.name}")
        print(f"X range     : {self.x_min} to {self.x_max}")
        print(f"Y range     : {self.y_min} to {self.y_max}")
        print(f"Materials   : {self.material}")
        print(f"Size        : {self.width()} x {self.height()} cm")
        print(f"Area        : {self.area()} cm^2")

# ============================================================
# PART 2: GeometryEngine Class
# ============================================================
class GeometryEngine:
    """
    Manages all 2D rectangular regions. 
    Validate geometry and builds the material mesh. 
    """
    def __init__(self, domain_x, domain_y):
        """_summary_

        Args:
            domain_x (float): total domain width in cm (Lx)
            domain_y (float): total domain width in cm (Ly)
        """
        self.domain_x   = domain_x      # store domain_x
        self.domain_y   = domain_y      # store domain_y 
        self.regions    = []            # empty list - regions added later
        
        # These get filled when build_mesh() is called
        self.material_map   = None      # 2D array of material names
        self.x_centers      = None      # 1D array of x cell centers 
        self.y_centers      = None      # 1D array of y cell centers 
        
    def add_region(self, region):
        """ 
        Add a Region object to the geometry.
        """
        self.regions.append(region)
        
    def get_region_at(self, x, y):
        """ 
        Return the region containing point (x,y)
        Returns None if no region found.
        """
        for region in self.regions:
            if region.contains(x,y):
                return region 
        return None 
    
    def build_mesh(self, nx, ny):
        """Build 2D material map using cell centers

        Args:
            nx (int): number of cells in x direction
            ny (int): number of cells in y direction
        """
        
        # Step 1: Calculate cell sizes 
        dx = self.domain_x /nx 
        dy = self.domain_y /ny 
        
        # Step 2: Build cell center array
        # x centers: dx/2, 3dx/2, 5dx/2...
        self.x_centers  = np.array([(i+0.5)*dx for i in range(nx)])
        self.y_centers  = np.array([(j+0.5)*dy for j in range(ny)])
        
        # Step 3: build empty material map 
        self.material_map   = np.empty((ny,nx), dtype=object)
        
        # Step 4: fill the material map 
        for j in range(ny):
            for i in range(nx):
                x       = self.x_centers[i]
                y       = self.y_centers[j]
                region  = self.get_region_at(x,y)
                
                if region is not None:
                    self.material_map[j][i] = region.material
                else: 
                    self.material_map[j][i] = "UNDEFINED"

        print(f" Mesh build: {nx} x {ny} cells")
        print(f" dx = {dx:.4f} cm, dy = {dy:.4f} cm")

    def summary(self):
        """Print the summary of the geometry"""
        print(f" Domain     : {self.domain_x} x {self.domain_y} cm")
        print(f" Regions    : {len(self.regions)}")
        for region in self.regions:
            region.describe()
            print()
        
# ============================================================
# PART 3: GeometryVisualizer Class
# ============================================================
class GeometryVisualizer:
    """ 
    Visualizes 2D reactor geometry with color-coded materials,
    dimension labels, mesh overlay, and interactive zoom.
    """
    
    MATERIAL_COLORS = {
        "water"       : "steelblue",
        "uo2_fuel"    : "orange",
        "mox_fuel"    : "tomato",
        "control_rod" : "dimgray",
        "void"        : "white",
        "UNDEFINED"   : "red",
    }
    
    def __init__(self, engine):
        """
        engine (GeometryEngine): the geometry to visualize
        """
        self.engine = engine 
    
    def _get_color(self,material_name):
        """Return color for a material, auto-assign if unknown"""
        if material_name in self.MATERIAL_COLORS:
            return self.MATERIAL_COLORS[material_name]
        
        # Auto assign a color for unkown materials
        auto_colors = list(mcolors.TABLEAU_COLORS.values())
        idx         = len(self.MATERIAL_COLORS)%len(auto_colors)
        return auto_colors[idx]
    
    def plot_geometry(self, show_mesh=False, show_labels=True,
                      figsize = (10,8)):
        """_summary_

        Args:
            show_mesh (bool, optional): overlay the numerical mesh grid
            show_labels (bool, optional): show region names and dimensions
            figsize (tuple, optional): figure size in inches
        """
        fig, ax = plt.subplots(1,1,figsize=figsize)

        drawn_materials = {}
        for region in self.engine.regions:
            color   = self._get_color(region.material)
            
            rect    = patches.Rectangle(
                (region.x_min, region.y_min),
                region.width(),
                region.height(),
                linewidth = 2,
                edgecolor = "black",
                facecolor = color,
                alpha = 0.7
            )
            ax.add_patch(rect)
            
            # track for legend (only add each material once)
            if region.material not in drawn_materials:
                drawn_materials[region.material] = color 
            
            # --- Step 2: add region label in center ---
            if show_labels:
                cx = (region.x_min+region.x_max)/2
                cy = (region.y_min+region.y_max)/2
                ax.text(cx,cy,
                        f"{region.name} \n ({region.material})",
                        ha = "center", va="center",
                        fontsize = 8, fontweight = "bold",
                        color = "white",
                        bbox = dict(boxstyle = "round,pad=0.2",
                                    facecolor="black", alpha = 0.4))
        # Step 3: Draw mesh grid overlay
        if show_mesh and self.engine.x_centers is not None:
            # Draw vertical lines at cell bounaries
            dx = self.engine.x_centers[1]-self.engine.x_centers[0]
            dy = self.engine.y_centers[1]-self.engine.y_centers[0]
            
            # x boundaries
            x_edges = [self.engine.x_centers[0]-dx/2]
            for xc in self.engine.x_centers:
                x_edges.append(xc+dx/2)
                
            # y boundries
            y_edges = [self.engine.y_centers[0]-dy/2]
            for yc in self.engine.y_centers:
                y_edges.append(yc+dy/2)
                
            for xe in x_edges:
                ax.axvline(x=xe, color="white",
                           linewidth=0.5, linestyle="--", alpha=0.6)
                
            for ye in y_edges:
                ax.axhline(y=ye, color="white",
                           linewidth=0.5, linestyle="--", alpha=0.6)
                
        # --- Step 4: Add dimension arrows on the edges ---
        if show_labels:
            for region in self.engine.regions:
                # Label width at the bottom of each region
                ax.annotate(
                    f"{region.width():.0f} cm",
                    xy = ((region.x_min+region.x_max)/2, region.y_min),
                    xytext = (0,-20),
                    textcoords = "offset points",
                    ha="center", fontsize = 7, color="black"
                )
        # --- Step 5: Build legend from draw materials ---
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=color, edgecolor="black",
                  label = mat, alpha=0.7)
            for mat, color in drawn_materials.items()
        ]
        
        ax.legend(handles=legend_elements, loc="upper right", fontsize=7)
        
        # --- Step 6: Final formatting ---
        ax.set_xlim(0, self.engine.domain_x)
        ax.set_ylim(0, self.engine.domain_y)
        ax.set_xlabel("x (cm)", fontsize=11)
        ax.set_ylabel("y (cm)", fontsize=11)
        ax.set_title("ANGKOR - 2D Geometry Map", fontsize=13)
        ax.set_aspect("equal")
        plt.tight_layout()
        plt.show()
            
    

# ============================================================
# TEST
# ============================================================
if __name__ == "__main__":
    print("=" * 50)
    print("  TESTING: Region Class")
    print("=" * 50)

    # Test 1: Create valid regions
    fuel = Region(
        name="fuel_center",
        x_min=20.0, x_max=80.0,
        y_min=20.0, y_max=80.0,
        material="uo2_fuel"
    )
    reflector = Region(
        name="reflector",
        x_min=0.0,  x_max=20.0,
        y_min=0.0,  y_max=100.0,
        material="water"
    )

    print("\n--- Fuel Region ---")
    fuel.describe()

    print("\n--- Reflector Region ---")
    reflector.describe()

    # Test 2: contains()
    print("\n--- Testing contains() ---")
    print(f"  fuel.contains(50, 50) : {fuel.contains(50, 50)}")   # True
    print(f"  fuel.contains(10, 10) : {fuel.contains(10, 10)}")   # False
    print(f"  fuel.contains(20, 20) : {fuel.contains(20, 20)}")   # True

    # Test 3: bad region triggers error
    print("\n--- Testing Validation ---")
    try:
        bad = Region("bad", 80.0, 20.0, 0.0, 100.0, "water")
    except ValueError as e:
        print(f"  Caught correctly: {e}")
    print("  Test finished")
    
    
    #-------------------------------------------------------
    print("\n" + "=" * 50)
    print("  TESTING: GeometryEngine")
    print("=" * 50)

    # Build a simple PWR-like 1D slice geometry
    engine = GeometryEngine(domain_x=100.0, domain_y=100.0)

    engine.add_region(Region("left_reflector",  0.0,  20.0,
                              0.0, 100.0, "water"))
    engine.add_region(Region("fuel",           20.0,  80.0,
                              0.0, 100.0, "uo2_fuel"))
    engine.add_region(Region("right_reflector",80.0, 100.0,
                              0.0, 100.0, "water"))

    print("\n--- Geometry Summary ---")
    engine.summary()

    print("\n--- Building Mesh ---")
    engine.build_mesh(nx=100, ny=100)

    print("\n--- Material Map ---")
    print(engine.material_map)

    print("\n--- Spot Checks ---")
    print(f"  Cell (0,0) material : {engine.material_map[0][0]}")  # water
    print(f"  Cell (5,5) material : {engine.material_map[5][5]}")  # uo2_fuel
    print(f"  Cell (9,9) material : {engine.material_map[9][9]}")  # water

    # ------Test visulization ----
    print("\n" + "=" * 50)
    print("  TESTING: GeometryVisualizer")
    print("=" * 50)

    viz = GeometryVisualizer(engine)

    print("\n--- Plot 1: Geometry only ---")
    viz.plot_geometry(show_mesh=False, show_labels=True)

    print("\n--- Plot 2: Geometry + Mesh overlay ---")
    viz.plot_geometry(show_mesh=True, show_labels=True)