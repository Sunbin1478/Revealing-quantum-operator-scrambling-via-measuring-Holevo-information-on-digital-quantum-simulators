# Quantum Operator Scrambling via Holevo Information Measurement

The main code and data for our research are available in the paper: [Revealing quantum operator scrambling via measuring Holevo information on digital quantum simulators](https://doi.org/10.48550/arXiv.2503.00918).

## Repository Structure

### Data Directories
- **results**: Contains all image data presented in the article
  - **Recording_of_subcircuit**: Corresponding to Figure 2
  - **noiseless_results**: Corresponding to Figures 3, 4, and 5
  - **mitigation_error_results**: Corresponding to Figure 6

### Environment
- **venv**: Complete environment configuration used by the code

### Trace Distance
- **Sample_random_samples** : 
  This program generates mixed states at different system sizes according to random subcircuits and records the sequences of quantum operational gates.
 - **plot_trace_distance**:
  This program reads the recorded data and plots the trace distance of the mixed state and the maximum mixed state with respect to the number of samples, M.
### Spatial-temporal Pattern
- **Pauli_operator_scrambling**：you can adjust the initial operator and conditions to get details of the spread of operators on the lattice points.
- **Calculate_holevo_information**：calculates Holevo information based on the evolutionary details of the initial operators you want to compare
- **heat_maps**:visualizing spatial-temporal pattern.
  

### Noise And Error Mitigation
We provide two implementation approaches:
- **HolevoinformationIX_noise_random**: Our proposed efficient approximation method
- **HolevoinformationXY_noise_Bell**: The exact measurement method
Both methods can be adjusted for specific experimental details as needed.



