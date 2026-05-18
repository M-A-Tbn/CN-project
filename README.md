# Complex Networks Project: Galaxy Cluster Identification

This repository contains the codebase and data analysis for a course project in **Complex Networks**, focusing on identifying and analyzing galaxy groups and clusters. 

## Project Overview

The project aims to construct a network representation of galaxies and apply complex network methodologies (such as percolation scanning and susceptibility analysis) to detect galaxy clusters. The results are then evaluated against the ground-truth labeled galaxy clusters from the Tempel et al. (2014) SDSS DR10 catalog.

### Key Phases:
1. **Data Acquisition:** Querying the Tempel et al. SDSS DR10 catalog via Vizier to obtain a volume-limited main galaxy sample.
2. **Network Construction & Percolation:** Applying linking lengths (both perpendicular and parallel) to build the network and compute topological observables like susceptibility to find the critical phase transition.
3. **Evaluation:** Comparing the detected galaxy groups with the Tempel reference catalog using metrics such as completeness and purity across different richness bins and redshift slices.

## Repository Structure

- `Notebook/CN-Project.ipynb`: The core Jupyter Notebook containing the data pipeline, network construction, and evaluation code.
- `Notebook/requirements.txt`: Python dependencies required to run the notebook.
- `Report/`: LaTeX source files and templates for the final project report.

## Getting Started

### Prerequisites

Ensure you have Python installed. It is recommended to use a virtual environment.

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/M-A-Tbn/CWN-project.git
   cd CWN-project
   ```

2. Install the required dependencies:
   ```bash
   pip install -r Notebook/requirements.txt
   ```

### Running the Code

Navigate to the `Notebook` directory and start Jupyter to explore the analysis:
```bash
cd Notebook
jupyter notebook CN-Project.ipynb
```

## Feedback and Collaboration

Feedback and discussions are highly encouraged! Please refer to the **Discussions** tab on GitHub if you are a professor, advisor, or interested researcher wanting to provide guidance on the methodology, code, or draft report.
