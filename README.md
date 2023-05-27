# Pairwise Interaction Energy Calculator

This script is designed to calculate the pairwise interaction energy between residues within a cutoff distance of 12 angstroms during molecular dynamics simulations. It provides efficient and fast computation of the van der Waals energy, electrostatic energy, and total energy for each frame of the simulation. 

## Features

- Calculates pairwise interaction energies between residues within a specified cutoff distance.
- Outputs a pickle file containing van der Waals energy, electrostatic energy, and total energy for each frame.
- Demonstrates a high correlation of 99% with the results obtained from the NAMD molecular dynamics engine.
- Provides significant speed improvements for energy calculations compared to other methods.

## Usage

1. Ensure that you have the required input files, including the protein structure file (PSF) and the molecular dynamics trajectory file (DCD).
2. Set the appropriate file paths for the PSF and DCD files in the script.
3. Run the script to perform the pairwise energy calculations.
4. The output will be saved as a pickle file, which can be easily loaded and analyzed using Python.
5. Analyze the results, including the van der Waals energy, electrostatic energy, and total energy at each frame.

## Requirements

- MDAnalysis library
- NumPy library
- Parmed library (CharmmPsfFile, CharmmParameterSet)

## Performance

This script has been optimized for speed, allowing for fast calculation of pairwise interaction energies. The efficient implementation results in significantly faster computations compared to other methods commonly used in the field.

## Note

Please note that the accuracy and performance of the script are dependent on the quality and reliability of the input files. Ensure that the PSF and DCD files are correctly formatted and represent the desired molecular system accurately.

Feel free to explore this script and utilize its capabilities for your molecular dynamics studies. If you encounter any issues or have suggestions for improvements, please let us know.
