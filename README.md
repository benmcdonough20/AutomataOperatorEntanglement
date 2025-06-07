# Code for: Bridging Classical and Quantum Information Scrambling with the Operator Entanglement Spectrum

This repository contains the code used to produce the results presented in our paper:

**Bridging Classical and Quantum Information Scrambling with the Operator Entanglement Spectrum**  
Ben T. McDonough et al. 
https://arxiv.org/abs/2505.05575

## Overview

The scripts and notebooks in this repository were used to generate the figures and data analysis in the main text and supplement.

## Repository Structure

```
.
├── README.md               # This file
├── run_all.sh              # Run all analysis and regenerate figures
├── environment/            # Julia environment
├── src/                    # Code for simulation and analysis 
│   ├── cluster/                # Numerical experiments run on the ISU compute cluster Pronto
│   ├── MPOs/                   # MPO simulations of local unitary and automaton circuits
│   ├── permutations/           # Numerical experiments for random automatons, automaton-evolved operators, and Bernoulli matrices
│   ├── evolution/              # Comparison of level spacings between classical matrix ensembles and evolved observables
│   └── produce_figures/        # Analysis code used to produce figures
├── data/                   # Output data (linked with GitHub LFS)
└── figures/                # Figures produced by the code
```

## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/benmcdonough20/AutomataOperatorEntanglement.git 
   cd AutomataOperatorEntanglement
   ```

2. Activate Julia environment and install dependencies:
   ```julia
   using Pkg; Pkg.activate("environment")
   Pkg.instantiate()
   ```

3. Run files in ```src/produce_figures``` to generate raw figures.

This will reproduce the core numerical results and figures from the manuscript.

## Dependencies

All of the dependencies are listed in the ```environment/``` folder, and can be installed by running ```Pkg.instantiate()```.

## Reproducibility

Some simulations may take several days on a standard laptop. Generated data is saved to the `data/` directory.

## License

This code is made available under the MIT License. See `LICENSE` for details.

## Contact

For questions or comments, please contact ```ben.mcdonough@colorado.edu``` or open an issue.

---

### Citation

If you use this code in your own work, please cite:

```latex
@article{mcdonough2025bridging,
  title={Bridging Classical and Quantum Information Scrambling with the Operator Entanglement Spectrum},
  author={McDonough, Ben T and Chamon, Claudio and Wilson, Justin H and Iadecola, Thomas},
  journal={arXiv preprint arXiv:2505.05575},
  year={2025}
}
```
