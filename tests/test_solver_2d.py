import sys, os 
sys.path.insert(0, os.path.join(
    os.path.dirname(__file__), '..', 'angkor'))

from input_reader   import InputReader
from solver_2d      import Solver2D

def test_solver_keff():
    """
    k-eff must match known result within tolerance.
    """
    reader  = InputReader("input/pwr_2d.yaml")
    reader.read()
    
    solver  = Solver2D(reader.engine, 
                       reader.materials,
                       reader.solver)
    solver.solve()
    assert abs(solver.k_eff - 1.179618)< 0.00010