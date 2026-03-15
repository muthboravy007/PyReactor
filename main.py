import sys 
import yaml 
import os 
import time 

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'angkor'))

from input_reader import InputReader
from solver_2d    import Solver2D
from output_2d    import Output2D

# ------------------------------------
# - BANNER 
# ------------------------------------

def print_banner():
    print("="*60)
    print("     ANGKOR ⚛️  🇰🇭")
    print("     Advanced Neutron Group Diffusion")
    print("     K-eigenvalue of Reactors")
    print() 
    print("     Institute of Technology of Cambodia (ITC)")
    print("     Open-Source reactor physics code")
    print("="*60)
    print()


def main(input_file):
    """ 
    Run a full ANGKOR Simulation from an input file. 
    Args: 
        input_file (str): path to YAML input file 
    """
    print_banner()
    
    # --- Check input file exist or not -------------------
    if not os.path.exists(input_file):
        print(f"    Error: Input file not found: {input_file}")
        print(f"    Usage: python main.py input/pwr_2d.yaml")
        sys.exit(1)
        
    print(f"    Input file : {input_file}")
    print() 
    
    # ---- Phase 1: Read Input ----------------------------
    print(" [1/3]   Reading input file...")
    t0 = time.time()
    reader = InputReader(input_file)
    reader.read() 
    
    t1 = time.time()
    print(f"    Done in {t1-t0:.2f}s")
    print() 
    
    # ---- Phase 2: Solve ---------------------------------
    print(" [2/3]   Running 2D diffusion solver...")
    t2 = time.time()
    
    solver = Solver2D(
        engine      = reader.engine,
        materials   = reader.materials,
        settings    = reader.solver
    )
    solver.solve()
    t3 = time.time()
    print(f"    Done in {t3-t2:.2f}s")
    
    # ---- Phase 3: Output --------------------------------
    print(" [3/3]   Generating output...")
    
    # output folder = same name as input file (without .yaml)
    base_name   = os.path.splitext(os.path.basename(input_file))[0]
    output_dir  = os.path.join("output", base_name)
    os.makedirs(output_dir, exist_ok = True)
    
    out = Output2D(solver, reader)
    out.print_report()
    out.plot_all(save_dir = output_dir, show = False)
    out.save_report(save_dir = output_dir)
    
    # ---- Done -------------------------------------------
    t_total = time.time()-t0 
    print()
    print("="*60)
    print(f"    ANGKOR completed in {t_total:.2f} seconds")
    print(f"    Results saved to: {output_dir}/")
    print(f"    k-eff = {solver.k_eff:.6f}")
    print("="*60)
    
if __name__ == "__main__":
    if len(sys.argv) <2:
        print_banner()
        print(" Error: No Input file provided!")
        print()
        print(" Usage:")
        print(" python main.py input/pwr_2d.yaml")
        sys.exit(1)
    
    input_file = sys.argv[1]
    main(input_file)