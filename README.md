# PyReactor

A deterministic neutron diffusion solver for nuclear reactor core analysis.

Built using the finite difference method with 2-group energy treatment 
and power iteration for k-eff calculation.

## Features

- 2-group neutron diffusion theory (fast + thermal)
- 1D slab geometry with multi-region support
- Power iteration eigenvalue solver
- YAML-based input file system
- Flux and power distribution plots
- **Validated: 0 pcm error vs. analytical 2-group solution**

## Quick Start

### 1. Install dependencies
pip install -r requirements.txt

### 2. Run the example reactor
python main.py input/reactor.yaml

### 3. Expected output
k-eff = 1.308664 (reflected PWR-like slab)

## Example Input File

```code
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

## Validation

| Benchmark | PyReactor | Reference | Error | Status |
|---|---|---|---|---|
| Bare slab (analytical) | 1.305447 | 1.305446 | 0 pcm | Agree |
| Moderated slab (MCNP 6.2) | 1.30866 | 1.32811 | 1,464 pcm | Expected |

Run benchmark:
python benchmarks/bare_slab/run_benchmark.py

## Project Structure

```code 
PyReactor/
├── pyreactor/
│   ├── materials.py     # Material cross-section library
│   ├── geometry.py      # 1D mesh builder
│   ├── solver_1d.py     # 2-group diffusion solver + power iteration
│   └── output.py        # Results plotting and reporting
├── benchmarks/          # Validation test cases
├── input/               # Example input files
├── main.py              # Entry point
└── requirements.txt
```

## Roadmap

- [x] 1D slab geometry
- [x] 2-group diffusion theory
- [x] Power iteration (k-eff)
- [x] Analytical benchmark validation
- [ ] OpenMC benchmark comparison
- [ ] 2D x-y geometry
- [ ] Burnup/depletion module
- [ ] Multi-group support

## Author

Built as an open-source alternative to commercial deterministic 
reactor physics codes such as CASMO/SIMULATE.

## License
MIT License
free to use, modify, and distribute.
