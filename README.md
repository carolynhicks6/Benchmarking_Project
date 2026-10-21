# BIOI 500 Benchmarking Project: Assessing accuracy, memory usage, and runtime of Autodock Vina, rDock, and GNINA for redocking and virtual screening
### Code overview
Autodock Vina, rDock, and GNINA are three popular tools for molecular docking. My code utilizes a snakemake file to prepare input files for each of these tools given a randomized subset of files from the CASF-2016 PDB benchmarking set. Once the snakemake file prepares inputs, Autodock and rDock runs are integrated into the workflow. Due to GPU constraints, a web-based colab script was used to run GNINA separate from the rest of the project workflow. 
### Initial project setup
First, clone the repo into your machine and navigate into the project directory. 
```
git clone https://github.com/carolynhicks6/Benchmarking_Project.git
cd Benchmarking_Project/
```
### Creating conda environment, installing packages, and running code
First, create and activate the conda environment. 
```
conda create --name docking_pipeline
conda activate docking_pipeline
```
When the environment is active, the terminal prompt should look like this:
```
(docking_pipeline) user/path_to_directory
```
Next, install the required packages. 
```
conda install -c conda-forge meeko vina rDock obabel mdanalysis pdb2pqr
```
In order to run the snakemake file, use the following command.
```
snakemake -s Snakefile.smk -c 1 --use-conda
```
To deactivate the conda environment, run the following command:
```
conda deactivate docking_pipeline
```
### Colab script: GNINA
A link to the original script can be found here:
https://colab.research.google.com/github/MolSSI-Education/iqb-2025/blob/main/04_Cheminfo_crash_course.ipynb 
Current edits to run code on 1bzc are shown here:
https://colab.research.google.com/drive/13wA-c4P7b5n2CX9TgKDV5iN-RPrlenSr#scrollTo=uRm7Jby4V6Pq

### Inputs
The initial input is shown in sample_data, with each subdirectory corresponding to a pdb code. The files in each subdirectory correspond to .pdb formats for the protein and binding pocket, .mol2 for the protein, ligand, and geometry optimized ligand, and .sdf for the ligand. The inputs for each section of analysis are shown in the snakemake. 

### Expected outputs
The outputs for Autodock Vina is a .pdbqt file, GNINA is a .sdf file, and rDock is a .sd file. These coordinate files contain ranked poses and scoring data. 
 
