# pyexocross
**pyexocross** is a Python wrapper to run the fortran code **exocross** and generate opacities in the **petitRADTRANS** format.

 **pyexocross** generates and organizes the required files for **exocross** and transforms its output to be used directly as  _line-by-line_ opacities in [petitRADTRANS](https://petitradtrans.readthedocs.io/en/latest/).

## Installation
To install **pyexocross** you can clone this repository and run the following command:
```bash
git clone https://github.com/DGonzalezPicos/pyexocross
cd pyexocross
pip install .
```
### Pre-requisites
- **exocross** must be downloaded and compiled in your system. You can find the installation instructions [here](https://exocross.readthedocs.io/en/latest/compile.html) or in the [pRT documentation](https://petitradtrans.readthedocs.io/en/latest/content/opa_add.html#preparing-exocross-opacities-for-petitradtrans).
- petitRADTRANS wavelength grid: [wlen_petitRADTRANS.dat](https://keeper.mpdl.mpg.de/f/357e92d4e0bb4aca9039/?dl=1)
    place this file in the `pyexocross` directory (next to `setup.py`)

## Usage
First create a directory to store the ExoMol files, in this example we will work with NaH:
```bash
mkdir -p NaH/input
```
Download these files (or equivalent for your molecule) from [ExoMol](https://www.exomol.com/data/molecules) and place them in the `NaH/input` directory:
- `23Na-1H__Rivlin.def`: File with metadata
- `23Na-1H__Rivlin.pf`: Partition function file
- `23Na-1H__Rivlin.states`: File with the energy levels
- `23Na-1H__Rivlin.trans`: File with the transitions

Now in a python script, import the ExoCross class and create an instance with the molecule of interest:
```python
from pyexocross.exocross import ExoCross
exo = ExoCross('NaH')
exo.set_output('NaH_rivlin') # name of linelist
```
Next select a temperature and pressure grid or create your own (examples in `PT_grid/`):
```python
exo.load_PT_grid(file='PT_grids/PTpaths_test.ls')
```
This basic grid contains two temperature and two pressure points, so it will generate 4 cross-section files. To generate the cross-sections run:
```python
path_to_excross = '/your/absolute/path'
exo.set_path_exocross(path_to_excross)
exo.xcross_grid(Nprocs=4) # Nprocs is the number of parallel processes
```
The output will be saved in the `NaH_rivlin` directory (as defined above). These cross-sections can be used directly with petitRADTRANS to generate synthetic spectra by simply placing the output folder in `input_data/opacities/lines/line_by_line`.

