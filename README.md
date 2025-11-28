# Source Code for: An empirical study of heuristics for an integrated bin packing and scheduling problem

This repository contains the source code for the heuristics for the 1BP-VPT.

## 1. Repository Content

The solver is written in C++ and implements the following algorithms: a **Constructive Procedure**, **Variable Neighborhood Descent (VND)**, **Large Neighborhood Search (LNS)**, and **Sequential Value Correction (SVC)**.

| File | Description |
| :--- | :--- |
| `Main.cpp` | **Main execution file** required to compile the solver. |
| `Io.h` | Functions for managing input/output of problem instances. |
| `Data.h` | Declaration of problem-specific data structures. |
| `Objective.h` | Calculation of the objective function. |
| `Rng.h` | Random number generator functions. |
| `Heu.h` | Implementation of the **Constructive Heuristic** procedure. |
| `Vnd.h` | Implementation of the **Variable Neighborhood Descent (VND)** metaheuristic. |
| `Lns.h` | Implementation of the **Large Neighborhood Search (LNS)** metaheuristic. |
| `Svc.h` | Implementation of the **Sequential Value Correction (SVC)** heuristic. |

**Accompanying Files:**

* `p.txt`: **Input parameters file** required for defining instance configuration.
* `Summary_02_08.xlsx`: Tabulated computational results over the benchmark instances.
---

## 2. Usage

### Dependencies
* The code is written in **C++** and requires a **C++20 compiler**

### Running the solver
* To run the solver on Windows: 
```bash
Solver.exe <instance_file> <parameter_file> <method>
```
Where: 
* `<instance_file>`: Path to a problem instance.
* `<parameter_file>`: Path to the `p.txt` file.
* `<method>`: {heu, lns, svc, svc_vnd}.

### Required file format: 
* `p.txt` requires the following format:
* w1: Weight associated with the number of bins.
* w2: Weight associated with Maximum tardiness. 

### Benchmark instances: 
* The istances are expected to be in this [format]([https://jump.dev/JuMP.jl/stable/](https://github.com/AndreaPizzuti/One-dimensional-bin-packing-with-pattern-dependent-processing-time/tree/main/INSTANCES)
* This is the standard adopted in Marinelli, [F., Pizzuti, A., Wu, W., Yagiura, M., 2025. One-dimensional bin packing with pattern-dependent processing time. European Journal of Operational Research 322, 3, 770–782.](https://www.sciencedirect.com/science/article/pii/S0377221724008920)
