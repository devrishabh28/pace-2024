# OCM Solver

This solver is a a student submission of PACE 2024 challenge and contains a C++ implementation for solving One-Sided Crossing Minimization Prolem. The OCM problem involves arranging the nodes of a bipartite graph on two layers, with one layer fixed, in order to minimize the number of edge crossings. This problem is fundamental in graph drawing and is known to be NP-hard.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Usage](#usage)
- [License](#license)

## Introduction

The `solver.cpp` program implements various algorithms for handling graph data structures and calculating metrics such as the crossing number. The main class, `PaceGraph`, includes methods for loading graphs from files, counting crossings, and applying different heuristic methods to optimize the graph layout.

## Features

- Segment Tree implementation for efficient range queries and updates.
- Graph class (`PaceGraph`) with methods to load from and save to `.gr` files.
- Heuristic methods for optimizing graph layouts:
  - Median heuristic
  - Barycenter heuristic
  - Local search methods

## Usage

To build the solver, use the following command:

```bash
g++ solver.cpp -o solver
```

To run the solver, use the following command:

```bash
./solver
```

After running the solver, the graph needs to be inputted via stdin. One such example of an input is:<br>
p ocr 5 6 10<br>
1 6<br>
1 9<br>
1 11<br>
2 6<br>
2 10<br>
3 10<br>
4 10<br>
4 7<br>
4 8<br>
5 8<br>

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
