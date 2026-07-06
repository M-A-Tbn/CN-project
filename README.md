# Complex networks project: galaxy group identification

This repository contains the codebase and data analysis for a course project in **Complex Networks** to identify and analyze galaxy groups. 

## Project overview

The codebase constructs a network representation of galaxies and applies complex network methodologies, including percolation scanning and susceptibility analysis, to detect galaxy groups. The results are evaluated against the ground-truth labeled galaxy groups from the Tempel et al. (2014) SDSS DR10 catalog.

### Key phases

1. **Data acquisition:** Querying the Tempel et al. SDSS DR10 catalog through Vizier to obtain a volume-limited main galaxy sample.
2. **Network construction and percolation:** Application of perpendicular and parallel linking lengths to construct the network and compute topological observables, such as susceptibility, to determine the critical phase transition.
3. **Evaluation:** Comparison of the detected galaxy groups with the Tempel reference catalog using completeness and purity metrics across richness bins and redshift slices.

## Repository structure

- `Notebook/CN-Project.ipynb`: The core Jupyter Notebook containing the data pipeline, network construction, and evaluation code.
- `Notebook/requirements.txt`: Python dependencies required to run the notebook.
- `Report/`: LaTeX source files and templates for the final project report.

## Getting started

### Prerequisites

Ensure Python is installed. We recommend using a virtual environment.

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/M-A-Tbn/CN-Project.git
   cd CN-Project
   ```

2. Install the required dependencies:
   ```bash
   pip install -r Notebook/requirements.txt
   ```

### Running the code

Navigate to the `Notebook` directory and start Jupyter to explore the analysis:
```bash
cd Notebook
jupyter notebook CN-Project.ipynb
```

## Feedback and collaboration

Feedback and discussions are welcome. Please use the GitHub **Discussions** tab to provide guidance on the methodology, code, or draft report.

