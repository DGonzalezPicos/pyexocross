import pathlib
import numpy as np
import logging
# set up logging
logging.basicConfig(level=logging.INFO)

import subprocess as sp

class ExoCross:
    
    # path = pathlib.Path('data2/dario/exocross') # path to exocross executable
    path = pathlib.Path('/home/dario/phd/pyexocross') # path to exocross executable 
    path_exocross = pathlib.Path('/home/dario/phd/pyexocross/exocross')
    # TODO: make this a environment variable
    
    logger = logging.getLogger('ExoCross')  # Create a class-level logger
    # logger.setLevel(logging.INFO)  # Set logger level to INFO

    def __init__(self, name, resolution=1e6):
        
        self.name = name
        self.resolution = resolution

        # input files should be in a folder with the same name as the molecule
        self.input = self.path / f'{self.name}/input'
        self.input.mkdir(parents=True, exist_ok=True)
        
        # create a temporary folder for intermediate files
        self.tmp = self.path / f'{self.name}/tmp'
        self.tmp.mkdir(parents=True, exist_ok=True)
        
        
        
    def debug(self):
        self.logger.setLevel(logging.DEBUG)
        return self
        
    def load_PT_grid(self, file='PTpaths.ls'):
        """ Load the pressure and temperature grid from a file used by pRT:
        3 columns = (P, T, cross sections files)"""
        self.P_grid, self.T_grid, _ = np.genfromtxt(file).T
        return self
    
    def check_input(self):
        assert hasattr(self, 'input'), 'No input folder found'
        files = sorted(self.input.glob('*'))
        self.logger.info(f' Found {len(files)} files in input folder')
        # self.logger.debug(f'Files: {files}')
        suffixes = [file_i.suffix for file_i in files]
        # chekck the following extensions exist
        required = ['.def', '.pf', '.states', '.trans']
        assert all([ext in suffixes for ext in required]), f'Not all required files found: {required}'

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
        
    def generate_inp(self,
                     temperature,
                     pressure, 
                     Nprocs=4):
        
        if not hasattr(self, 'files'):
            self.check_input()
            
        # pffile = self.pf.name
        # states = self.states.name
        # transitions = [t.name for t in self.trans]
        # use absolute paths instead
        pffile = str(self.pf)
        states = str(self.states)
        transitions = [str(t) for t in self.trans]
        
        # out_file = f"{self.name}_{temperature:f}K_{pressure:f}bar.out"

        species_str = '\n 0 gamma 0.06 n 0.5 t0 296 ratio 1.\nend'

        # write the .inp file
        inp_file = self.input / f'{self.name}_input_{temperature:.1f}_{pressure:.0e}.inp'
        out_file = self.tmp / inp_file.name.replace('.inp', '.out')
        
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
            # close the file
            f.close()
        self.logger.info(f' Wrote {inp_file}')
        
        return inp_file, out_file      
    
    
    def generate_PTpaths(self, files_path, overwrite=False):
        
        if not isinstance(files_path, pathlib.Path):
            files_path = pathlib.Path(files_path)
            
        PTpaths = files_path / 'PTpaths.ls'

        assert files_path.exists(), f'Path {files_path} does not exist'
        
        if PTpaths.exists() and not overwrite:
            self.logger.info(f'{PTpaths} already exists. Use overwrite=True to overwrite')
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
        self.logger.info(f'Wrote {PTpaths}')
        return self
    
    
    def xcross(self, temperature, pressure, Nprocs=4):
        
        # create input file
        inp_file, out_file = self.generate_inp(temperature, pressure, Nprocs)
        print(out_file)
        # run exocross
        input_file = str(self.input / inp_file)
        output_file = str(self.tmp / out_file)
        self.logger.info(f'Running exocross with input file: {input_file}')
        self.logger.info(f'Output will be saved to: {output_file}')
        # command = [f"{self.path_exocross}/xcross.exe",
        #            f"< {input_file}",
        #            f"> {output_file}"]
        command = [f"{self.path_exocross}/xcross.exe"]

        arguments = ["<", input_file, ">", output_file]

        # call the command and send output to the logger
        # print(command+arguments)
        # self.logger.info(f'Running: {command}')
        # pass input and output to the command
        # sp.run(command+arguments)
        result = sp.run(["sh", "-c", f"{command[0]} < {input_file} > {output_file}"], capture_output=True)
        if result.returncode != 0:
            print("Error running command!")
            print(result.stderr.decode("utf-8"))  # Decode and print error messages
        # result = sp.run(command, capture_output=True, stdin=sp.PIPE)
        # # Check for successful execution (return code 0)
        # if result.returncode == 0:
        #     return result.returncode, result.stdout.decode()
        # else:
        #     raise RuntimeError(f"xcross.exe exited with non-zero code: {result.returncode}")
        
        print('done')
        self.logger.info(f'Finished running {command}')
        # return the output file
        
        return self
    
    

        
if __name__ == '__main__':
    
    exo = ExoCross('NaH').debug()
    exo.check_input()
    exo.read_definitions_file()
    # exo.generate_PTpaths('/home/dario/phd/pRT_input/input_data/opacities/lines/line_by_line/CO2_main_iso')
    # exo.load_PT_grid()
    # inp_file, out_file =exo.generate_inp(temperature=300, pressure=1e-5)
    exo.xcross(temperature=300, pressure=1e-5, Nprocs=1)
    
    