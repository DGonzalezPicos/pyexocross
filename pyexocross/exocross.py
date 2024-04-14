import pathlib
import os
import numpy as np
import logging
# set up logging
logging.basicConfig(level=logging.INFO)

import subprocess as sp
# import pool for multiprocessing
from multiprocessing import Pool

class ExoCross:
    
    # path = pathlib.Path('data2/dario/exocross') # path to exocross executable
    # path = pathlib.Path('/data2/dario/pyexocross/') 
    path_exocross = pathlib.Path('/data2/dario/pyexocross/exocross') # path to exocross executable 
    # TODO: make this a environment variable
    
    logger = logging.getLogger('ExoCross')  # Create a class-level logger
    # logger.setLevel(logging.INFO)  # Set logger level to INFO

    def __init__(self, name, resolution=1e6, path=None):
        
        self.name = name
        self.resolution = resolution
        
        path = path if path is not None else os.getcwd()
        self.path = pathlib.Path(path)
        self.logger.debug(f' Working in {self.path}')

        # input files should be in a folder with the same name as the molecule
        self.input = self.path / f'{self.name}/input'
        self.input.mkdir(parents=True, exist_ok=True)
        
        # create a temporary folder for intermediate files
        self.tmp = self.path / f'{self.name}/tmp'
        self.tmp.mkdir(parents=True, exist_ok=True)
        
        self.output = self.tmp # default output folder
        
        
        # check if the input folder exists and contains the required files
        self.check_input()
        self.read_definitions_file()
        
    def debug(self):
        self.logger.setLevel(logging.DEBUG)
        return self
    def set_path_exocross(self, path):
        self.path_exocross = pathlib.Path(path)
        # save as hidden file in the working directory
        with open('.path_exocross', 'w') as f:
            f.write(str(self.path_exocross))
        return self
    def find_path_exocross(self):
        ''' Find the path to the exocross executable'''
        if not hasattr(self, 'path_exocross'):
            if '.path_exocross' in os.listdir():
                with open('.path_exocross', 'r') as f:
                    self.path_exocross = pathlib.Path(f.read())
            else:
                self.logger.warning(' No exocross path found. Set it with `set_path_exocross()`')
        return self
        
    def load_PT_grid(self, file='PTpaths.ls'):
        """ Load the pressure and temperature grid from a file used by pRT:
        3 columns = (P, T, cross sections files)"""
        self.P_grid, self.T_grid, _ = np.genfromtxt(file).T    
        return self
    
    def check_input(self):
        assert hasattr(self, 'input'), 'No input folder found'
        files = sorted(self.input.glob('*'))
        assert len(files) > 0, f'No files found in {self.input}'
        # self.logger.info(f' Found {len(files)} files in input folder')
        # self.logger.debug(f'Files: {files}')
        suffixes = [file_i.suffix for file_i in files]
        # chekck the following extensions exist
        required = ['.def', '.pf', '.states', '.trans']
        assert all([ext in suffixes for ext in required]), f'Not all required files found: {required}'
        self.logger.info(f' Found all required files in {self.input}')
        self.files = files
        # find the files of the required extensions
        for ext in required:
            attr = [file_i for file_i in files if file_i.suffix == ext]
            if ext == '.trans': # allow for several transitions files
                setattr(self, ext[1:], attr)
            if ext == '.def':
                setattr(self, 'defi', attr[0]) # `def`` is not an allowed attribute name
            else:
                setattr(self, ext[1:], attr[0])
                
        self.trans = np.atleast_1d(self.trans)
        return self
    
    def read_definitions_file(self):
        assert hasattr(self, 'defi'), 'No definition file found'
        with open(self.defi, 'r') as f:
            lines = f.readlines()
            # separate the lines in the file by the '#' character
            sections = [line.split('#') for line in lines]
            # remove commas
            sections = [[s.strip() for s in section] for section in sections]
            # create a dictionary with the sections
            self.info = {str(sec[1]): sec[0] for sec in sections}
            
        self.label = f"{self.info['Iso-slug']}__{self.info['Isotopologue dataset name']}"
        self.mass = float(self.info['Isotopologue mass (Da) and (kg)'].split(' ')[0])
        return self        
    
    def set_broadening_params(self):
        # TODO: implement this
        pass
    
    def set_output(self, label=None):
        
        label = self.label.lower() if label is None else label
        self.output = self.path / f'{self.name}/{label}'
        self.output.mkdir(parents=True, exist_ok=True)
        
        self._create_molparam_id() # default molparam_id file
        return self
        
    def generate_inp(self,
                     temperature,
                     pressure, 
                     Nprocs=4):
        
        if not hasattr(self, 'files'):
            self.check_input()
            
        # use absolute paths 
        pffile = str(self.pf)
        states = str(self.states)
        transitions = [str(t) for t in self.trans]
        
        # write the .inp file
        inp_file = self.input / f'{self.name}_input_{temperature:.1f}_{pressure:.0e}.inp'
        out_file = self.tmp / inp_file.name.replace('.inp', '.out')
        
        # default broadening
        # TODO: enable user to change this in 'set_broadening_params()'
        species_str = '\n 0 gamma 0.06 n 0.5 t0 296 ratio 1.\nend'

        with open(inp_file, 'w') as f:
            f.write(f'Nprocs {Nprocs}\n')
            f.write('absorption\n')
            f.write('voigt\n')
            f.write('verbose 3\n')
            f.write('offset 60.\n')
            f.write(f'mass {int(np.round(self.mass,0))}\n')
            f.write(f'temperature {temperature:.1f}\n')
            f.write(f'pressure {pressure:.0e}\n')
            f.write('range 39. 91000.\n')
            f.write(f'R {self.resolution:.0f}\n')
            f.write(f'pffile {pffile}\n')
            f.write(f'output {out_file}\n')
            f.write(f'states {states}\n')
            f.write(f'transitions\n')
            for trans in transitions:
                f.write(f'\t"{trans}"\n')
                    
            f.write('end\n')
            f.write(f'species {species_str}')
            f.close()
        self.logger.debug(f' Wrote {inp_file}')
        
        return inp_file, out_file      
    
    
    def generate_PTpaths(self, files_path=None, overwrite=False):
        """ Generate a file with the pressure and temperature grid for the cross sections
        from the file names in `files_path`"""
        
        files_path = self.output if files_path is None else files_path
        if not isinstance(files_path, pathlib.Path):
            files_path = pathlib.Path(files_path)
            
        PTpaths = files_path / 'PTpaths.ls'

        assert files_path.exists(), f'Path {files_path} does not exist'
        
        if PTpaths.exists() and not overwrite:
            self.logger.info(f' {PTpaths} already exists. Use overwrite=True to overwrite')
            return self
            
        
        self.sigmas = sorted(files_path.glob('sigma*'))
        assert len(self.sigmas) > 0, 'No sigma files found'
        # P, T = [], []
        self.rows = []
        for f in self.sigmas:
            # self.logger.debug(f)
            t = f.name.split('.K')[0].split('_')[-1]
            p = f.name.split('.K')[1].split('bar')[0][1:]
            self.rows.append([p, t, f.name])
            # print
        p_float = [float(row[0]) for row in self.rows]
        t_float = [float(row[1]) for row in self.rows]
        # sort by pressure and temperature
        idx = np.lexsort((t_float, p_float))
        self.P_grid = np.array(p_float)[idx]
        self.T_grid = np.array(t_float)[idx]
        self.rows = np.array(self.rows)[idx]
        # save the PT paths to a file with extension 6 digits for the first two columns
        np.savetxt(PTpaths, self.rows, fmt=['%10s', '%10s', '%s'])
        self.logger.info(f' Wrote {PTpaths}')
        return self
    
    
    def xcross(self, temperature, pressure, Nprocs=1):
        """ Main function to run exocross with the given temperature and pressure"""
        
        # create input file
        inp_file, out_file = self.generate_inp(temperature, pressure, Nprocs)
        input_file = str(self.input / inp_file)
        output_file = str(self.tmp / out_file)
        self.logger.info(f' [xcross.exe] {input_file}')
        
        # find the path to the exocross executable
        self.find_path_exocross()
        command = f"{self.path_exocross}/xcross.exe"

        # call the command
        result = sp.run(["sh", "-c", f"{command} < {input_file} > {output_file}"], capture_output=True)
        if result.returncode != 0:
            print(result.stderr.decode("utf-8"))  # Decode and print error messages
        
        self.logger.info(f' --> finished {input_file}')        
        return self
    
    def xcross_grid(self, Nprocs=4):
        """ Run exocross for all the temperatures and pressures in the grid"""
        
        assert hasattr(self, 'T_grid'), 'No temperature grid found'
        # print section header
        print(f' \n *** xcross grid ***')
        
        if Nprocs > 1:
            # Create a pool of worker processes
            with Pool(processes=Nprocs) as pool:
                # Use pool.starmap to unpack arguments for each run_xcross call
                pool.starmap(self.xcross, zip(self.T_grid, self.P_grid))
            
        else:
            for T, P in zip(self.T_grid, self.P_grid):
                self.xcross(T, P)
        
        self.rebin_to_pRT()
        self.make_short()
        self.generate_PTpaths()
        return self
    
    def rebin_to_pRT(self):
        """ Rebin the cross sections to the pRT grid
        
        Adapted from pRT documentation"""
        
        files = sorted(self.tmp.glob('*.out.xsec'))
        self.logger.info(f' Found {len(files)} files in tmp folder')
        
        wave_pRT = np.genfromtxt('wlen_petitRADTRANS.dat')
        for i, file_i in enumerate(files):
            self.logger.info(f' Rebinning {file_i} ({i}/{len(files)})')
            T_str, P_str = file_i.stem.split('.out')[0].split('_')[-2:]
            
            T = float(T_str)
            P = float(P_str)
            
            # read the cross sections
            dat = np.genfromtxt(file_i)
            wavelength = 1./dat[:,0]
            sigma = dat[:,1]

            # Invert them to go from a accending wavenumber ordering
            # to an accending wavelength ordering.
            wavelength = wavelength[::-1]
            sigma = sigma[::-1]
            
            # interpolate to the pRT grid
            sig_interpolated_petit = np.interp(wave_pRT, wavelength, sigma)
                        
            new_file_i = self.tmp / f'sigma_{T:.0f}.K_{P:.6f}bar.dat'
            
            # save rebinned cross sections
            np.savetxt(new_file_i, np.column_stack((wave_pRT, sig_interpolated_petit)))
            self.logger.info(f' Wrote {new_file_i}')
        return self
    
    def _create_short_stream_lambs_mass(self, overwrite=True):
        """Create required file for make_short.f90 to shorten the wavelength grid"""
        
        # file with three rows
        file_name = self.tmp / "short_stream_lambs_mass.dat"
        # check if the file already exists
        if file_name.exists() and not overwrite:
            self.logger.debug(f' {file_name} already exists')
            return self
        
        min_wave_cm = "0.3d-4"
        max_wave_cm = "28d-4"
        molecular_mass_amu = f'{np.round(self.mass, 0):.0f}d0'
        # write file
        self.logger.debug(f' Writing {file_name}')
        with open(file_name, 'w') as f:
            f.write("# Minimum wavelength in cm\n")
            f.write(f"{min_wave_cm}\n")
            f.write("# Maximum wavelength in cm\n")
            f.write(f"{max_wave_cm}\n")
            f.write("# Molecular mass in amu\n")
            f.write(f"{molecular_mass_amu}")
            f.close()
        return self
    
    def _create_molparam_id(self):
        ''' Create the molparam_id file (required by pRT)'''
        file_name = self.output / "molparam_id.txt"
        if file_name.exists():
            self.logger.debug(f' {file_name} already exists')
            return self
        with open(file_name, 'w') as f:
            f.write(f"#### Species ID (A2) format\n")
            f.write("06\n")
            f.write("#### molparam value\n")
            f.write("1.0")
            
        self.logger.debug(f' Wrote {file_name}')
        return self
            
            
    def make_short(self):
        ''' Shorten the wavelength grid to the pRT grid using the fortran code make_short.f90
        '''
        # 1) create the list of files
        files = sorted(self.tmp.glob('sigma*bar.dat'))
        files = [f.name for f in files] # keep only the name
        
        self.logger.info(f' Found {len(files)} files in output folder')
        np.savetxt(self.tmp / 'sigma_list.ls', files, fmt='%s')
        
        # 2) create short_stream_lambs_mass.dat
        self._create_short_stream_lambs_mass()
        
        # 3) create folder for the output
        sp.run(["mkdir", "short_stream"], capture_output=True, cwd=self.tmp)
        
        # 3) compile and run the fortran code
        # copy the make_short.f90 file to the output folder
        sp.run(["cp", f"{self.path_exocross}/make_short.f90", f"{self.tmp}"])
        # compile the fortran code
        sp.run(["gfortran", "-o", "make_short", "./make_short.f90"], cwd=self.tmp)
        # call the executable
        sp.run(["./make_short"], cwd=self.tmp)
        
        # copy output files to the output folder
        if self.output != self.tmp:
            sp.run(["sh", "-c", f"cp -r {self.tmp}/short_stream/* {self.output}"])
        
        return self

        
if __name__ == '__main__':
    
    exo = ExoCross('CN').debug()
    exo.check_input()
    exo.read_definitions_file()
    exo.set_output('CN_test_multi')

    exo.load_PT_grid(file='PT_grids/PTpaths_test.ls')
    exo.xcross_grid(Nprocs=8)

    
    