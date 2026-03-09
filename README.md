# ANGKOR ⚛️ 🇰🇭
### Advanced Neutron Group-diffusion K-eigenvalue Of Reactors

A deterministic neutron diffusion code system for nuclear reactor core analysis. Built using the finite difference method with multi-group energy treatment and power iteration for k-eff calculation.

> *Named after Angkor Wat — the greatest symbol of Khmer civilization —  
> because Cambodia's scientific legacy continues.*

**Developed in Cambodia 🇰🇭 as an open-source alternative to commercial deterministic codes such as CASMO/SIMULATE.**

---

## Features

- 2-group neutron diffusion theory (fast + thermal groups)
- 1D slab geometry with multi-region support
- **2D rectangular geometry engine with material map**
- **Interactive geometry visualizer (color-coded, zoom/pan)**
- Power iteration eigenvalue solver (k-eff)
- YAML-based human-readable input file system
- Flux and power distribution plots
- **Validated: 0 pcm error vs. analytical 2-group solution**

---

## Quick Start

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Run the example reactor
```bash
python main.py input/reactor.yaml
```

### 3. Expected output
```
k-eff = 1.308664  (reflected PWR-like slab)
```

### 4. Visualize 2D geometry
```bash
python angkor/geometry_2d.py
```

---

## Example Input File

```yaml
geometry:
  type: 1D_slab
  length_cm: 200.0
  mesh_points: 200

regions:
  - name: left_reflector
    start: 0.0
    end: 20.0
    material: water
  - name: fuel
    start: 20.0
    end: 180.0
    material: uo2_fuel
  - name: right_reflector
    start: 180.0
    end: 200.0
    material: water

materials:
  uo2_fuel:
    D1: 1.40
    D2: 0.37
    sigma_a1: 0.0095
    sigma_a2: 0.0820
    sigma_s12: 0.0210
    nu_sigma_f1: 0.0060
    nu_sigma_f2: 0.1350
  water:
    D1: 1.13
    D2: 0.16
    sigma_a1: 0.0003
    sigma_a2: 0.0210
    sigma_s12: 0.0460
    nu_sigma_f1: 0.0
    nu_sigma_f2: 0.0
```

---

## Validation

| Benchmark | ANGKOR | Reference | Error | Status |
|---|---|---|---|---|
| Bare slab (analytical) | 1.305447 | 1.305446 | 0 pcm | ✅ Verified |
| Moderated slab (MCNP 6.2) | 1.30866 | 1.32811 | 1,464 pcm | ✅ Expected* |

> *Diffusion theory vs Monte Carlo difference is physically expected.  
> Full transport correction planned in future versions.

Run benchmark:
```bash
python benchmarks/bare_slab/run_benchmark.py
```

---

## Project Structure

```
ANGKOR/
├── angkor/
│   ├── materials.py       # Material cross-section library
│   ├── geometry.py        # 1D mesh builder
│   ├── geometry_2d.py     # 2D rectangular geometry engine + visualizer
│   ├── solver_1d.py       # 2-group diffusion solver + power iteration
│   └── output.py          # Results plotting and reporting
├── benchmarks/            # Validation test cases
│   ├── bare_slab/
│   └── mcnp_benchmark/
├── docs/                  # Documentation
├── input/                 # Example input files
├── tests/                 # Unit tests
│   ├── test_geometry.py
│   ├── test_materials.py
│   └── test_solver_1d.py
├── main.py                # Entry point
├── requirements.txt
└── README.md
```

---

## Roadmap

### Version 1.0 — 1D Foundation 
- [x] 1D slab geometry
- [x] 2-group neutron diffusion theory
- [x] Power iteration (k-eff)
- [x] Analytical benchmark validation
- [x] YAML input file system

### Version 1.1 — 2D Geometry 
- [x] 2D rectangular geometry engine
- [x] Interactive geometry visualizer
- [x] Material color map with labels and mesh overlay

### Version 1.2 — 2D Solver (In Progress)
- [ ] 2D finite difference diffusion solver
- [ ] 2D flux and power distribution
- [ ] MCNP benchmark comparison (2D)

### Version 2.0 — Advanced Physics (Planned)
- [ ] Multi-group energy support (4-group, 8-group)
- [ ] Cylindrical geometry (fuel pins, control rods)
- [ ] Full core PWR / BWR / SMR models

### Version 3.0 — Burnup & Depletion (Future)
- [ ] Burnup calculation
- [ ] Isotope depletion (Bateman equations)
- [ ] Cycle length prediction

---

## About

ANGKOR is developed at the **Institute of Technology of Cambodia (ITC)** as part of a research initiative to build open-source nuclear reactor analysis tools accessible to the global nuclear engineering community.

---

## Author

**MUTH Boravy**  
Institute of Technology of Cambodia (ITC)  
Phnom Penh, Cambodia 🇰🇭

---

## License

MIT License — free to use, modify, and distribute.  
See [LICENSE](LICENSE) for details.
