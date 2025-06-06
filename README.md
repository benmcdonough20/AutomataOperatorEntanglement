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
├── requirements.txt        # Python dependencies
├── src/                    # Main simulation and analysis code
│   ├── circuits.py
│   ├── entanglement.py
│   └── plot_helpers.py
├── data/                   # Output data (optional or .gitignored)
├── figures/                # Reproduced figures from the paper
└── run_all.sh              # Script to regenerate all main figures
```

## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/operator-entanglement-spectra.git
   cd operator-entanglement-spectra
   ```

2. Create a virtual environment and install dependencies:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

3. Run the main analysis:
   ```bash
   bash run_all.sh
   ```

This will reproduce the core numerical results and figures from the manuscript.

## Dependencies

The code was developed using Python 3.10 and relies on:

- NumPy
- SciPy
- Matplotlib
- QuTiP (for some circuit simulations)
- JAX (optional: for faster simulation backends)

All dependencies are listed in `requirements.txt`.

## Reproducibility

The random seeds used for all stochastic sampling are fixed in the code. Some simulations may take several hours on a standard laptop. Intermediate data is saved to the `data/` directory.

## License

This code is made available under the MIT License. See `LICENSE` for details.

## Contact

For questions or comments, please contact [your email] or open an issue.

---

### Citation

If you use this code in your own work, please cite:

> Ben McD. et al., *Bridging Classical and Quantum Information Scrambling with the Operator Entanglement Spectrum*, submitted to *Physical Review X* (2025).
