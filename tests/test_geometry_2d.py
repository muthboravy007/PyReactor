import sys, os 
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..' ,'angkor'))

from geometry_2d import Region, GeometryEngine 

# ---- Test 1: Test region, whether we create correctly or not -----------
def test_region_properties():
    r = Region("fuel", 20.0, 80.0, 0.0, 100.0, "uo2_fuel")
    assert r.name       == "fuel"
    assert r.material   == "uo2_fuel"
    assert r.width()    == 60
    assert r.height()   == 100
    assert r.area()     == 6000
    
# ---- Test 2: Test contain() working or not 
def test_region_contains():
    r = Region("fuel", 20.0, 80.0, 0.0, 100.0, "uo2_fuel")
    
    assert r.contains(50,50) == True
    assert r.contains(10,50) == False
    assert r.contains(20,50) == True 
    
# =============================================
# TEST 3: bad region raises error
# =============================================
def test_region_invalid():
    import pytest
    with pytest.raises(ValueError):
        Region("bad", 80.0, 20.0,      
               0.0, 100.0, "water")

# =============================================
# TEST 4: GeometryEngine builds correctly
# =============================================
def test_engine_material_map():
    engine = GeometryEngine(100.0, 100.0)
    engine.add_region(Region("left",  0.0,  20.0, 0.0, 100.0, "water"))
    engine.add_region(Region("fuel", 20.0,  80.0, 0.0, 100.0, "uo2_fuel"))
    engine.add_region(Region("right",80.0, 100.0, 0.0, 100.0, "water"))
    engine.build_mesh(10, 10)

    # Cell (0,0) is at x=5cm → should be water
    assert engine.material_map[0][0] == "water"

    # Cell (5,5) is at x=55cm → should be fuel
    assert engine.material_map[5][5] == "uo2_fuel"
    