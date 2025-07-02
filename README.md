## Code to generate data, figures, and tables for "Fast Mølmer-Sørensen gates in trapped-ion quantum processors with compensated carrier transition" (https://arxiv.org/abs/2501.02387)

Below, there is a brief description of the code and data present in the repository.

* In the `figures` directory, the scripts `figure_n.py`, `n = 1..5` produce the figures numbered as in the manuscript.
* In the `tables` directory, the scripts `table_n.py`, `n = 1, 2` produce the $\LaTeX$-formatted tables numbered as in the manuscript.
* Both `figure_n.py` and `table_n.py` use the data files `ms_data.db` and `n_ions_dep.txt` and `.npz` archives  stored in the `./data` directory. Also, they use the code imported from the `timslib` Python package (see below).
* The `timslib` package source code is located in `src/timslib`, and the package can be installed from source in a virtual environment via `pip install`. It contains the code for analytical calculations of the quantum dynamics of trapped-ion chains illuminated by bichromatic amplitude-shaped laser pulses.
* The scripts to generate data files are in the `generate_data` directory. They contain functions for the numerical solution of the time-dependent Schrödinger equation for trapped ions, and they use the code from `timslib` for amplitude pulse shaping and related calculations.
