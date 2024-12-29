# LM-Merger

LM-Merger is a Python package designed for logic model merging. It provides utilities for input/output operations, parsing models, and merging them using three approaches: 'AND', 'OR', and 'Inhibitor Wins'.

## Requirements

This package is intended to be used within the **CoLoMoTo Interactive Notebook** Docker environment.

## Installation and Setup

### Step 1: Install the CoLoMoTo Docker Image
Ensure Docker is installed on your system. Install the required CoLoMoTo Docker image:

```bash
pip install -U colomoto-docker    # or python3 -m pip install ..
```
The CoLoMoTo notebook can be started by executing in a terminal:

```bash
colomoto-docker
```
Please refer to the [CoLoMoTo documentation](https://colomoto.github.io/colomoto-docker) for more information.

### Step 2: Install LM-Merger
0. Open a terminal in the CoLoMoTo Docker container.

1. Clone the LM-Merger repository into the Docker container:
```bash
git clone https://github.com/lunaticstarr/lmmerger.git
```

2. Install the required dependency (python-libsbml) and the LM-Merger package:
```bash
pip install -r requirements.txt
```

## Usage

To use LM-Merger, you can import and use the package in your Python script or Jupyter Notebook:

```python
from lmmerger import io, merger

## Load the networks
network1 = io.read_network(sbml_file1) # load the SBML-qual model
network2 = io.read_network(txt_file2) # can also be a txt or csv file
print("Network #1:")
print(network1)
print("Network #2:")
print(network2)

# Merge using OR approach
merged_network_or = merger.merge_networks([network1, network2], method="OR")
# Merge using AND approach
merged_network_and = merger.merge_networks([network1, network2], method="AND")
# Merge using Inhibitor Wins approach
merged_network_inhibitor = merger.merge_networks([network1, network2], method="Inhibitor Wins")

# Write the merged networks to files
merged_and_name = "merged_and_" + model1name + "_" + model2name
io.write_network_to_file(merged_network_and, home + "LogicModelMerger/Models/" + merged_and_name, format="text") # write to txt file
io.write_network_to_file(merged_network_and, home + "LogicModelMerger/Models/" + merged_and_name, format="sbml") # write to SBML-qual file

```

### Input and output formats
- Text file similar to the R package `Boolnet`, as described [here](https://rdrr.io/cran/BoolNet/man/loadNetwork.html)
- SBML-qual format
