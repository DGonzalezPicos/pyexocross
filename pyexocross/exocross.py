import pathlib
import numpy as np

from mass import table


class ExoCross:
    
    # path = pathlib.Path('data2/dario/exocross') # path to exocross executable
    path = pathlib.Path('/home/dario/phd/pyexocross') # path to exocross executable 
    # TODO: make this a environment variable
    
    def __init__(self, name, mass=None, resolution=1e6):
        
        self.name = name
        if mass is None:
            self.mass = table[name]
            print(f' Mass of {name} = {self.mass}')
            
        self.resolution = resolution
        
        
        self.input = self.path / f'{self.name}/input'
        self.input.mkdir(parents=True, exist_ok=True)
        # 
        
    def load_PT_grid(self, file='PTpaths.ls'):
        """ Load the pressure and temperature grid from a file used by pRT:
        3 columns = (P, T, cross sections files)"""
        self.P_grid, self.T_grid, _ = np.genfromtxt(file).T
        return self
    
    def check_input(self):
        assert hasattr(self, 'input'), 'No input folder found'
        files = sorted(self.input.glob('*'))
        print(f'Found {len(files)} files in input folder')
        print(files)
        suffixes = [file_i.suffix for file_i in files]
        # chekck the following extensions exist
        required = ['.pf', '.states', '.trans']
        assert all([ext in suffixes for ext in required]), f'Not all required files found: {required}'

        self.files = files
        # find the files of the required extensions
        for ext in required:
            attr = [file_i for file_i in files if file_i.suffix == ext]
            if ext == '.trans': # allow for several transitions files
                setattr(self, ext[1:], attr)
            else:
                setattr(self, ext[1:], attr[0])
        return self
    
        
    def generate_inp(self,
                     temperature,
                     pressure, 
                     Nprocs=4):
        
        if not hasattr(self, 'files'):
            self.check_input()
            
        # pffile = sorted(folder.glob('*.pf'))[0].name
        pffile = self.pf.name
        states = self.states.name
        transitions = [t.name for t in self.trans]
        
        output = f"{self.name}_{temperature:f}K_{pressure:f}bar.out"

        species_str = '\n 0 gamma 0.06 n 0.5 t0 296 ratio 1.\nend'

        # write the .inp file
        inp_out = self.input / f'{self.name}_input.inp'
        with open(inp_out, 'w') as f:
            f.write(f'Nprocs {Nprocs}\n')
            f.write('absorption\n')
            f.write('voigt\n')
            f.write('verbose 3\n')
            f.write('offset 60.\n')
            f.write(f'mass {int(self.mass)}\n')
            f.write(f'temperature {temperature}\n')
            f.write(f'pressure {pressure:f}\n')
            f.write('range 39. 91000.\n')
            f.write(f'R {self.resolution}\n')
            f.write(f'pffile {pffile}\n')
            f.write(f'output {output}\n')
            f.write(f'states {states}\n')
            f.write(f'transitions\n')
            for trans in transitions:
                f.write(f'\t"{trans}"\n')
                    
            f.write('end\n')
            f.write(f'species {species_str}')
            # close the file
            f.close()
        print(f'Wrote {inp_out}')
        
        return inp_out, output      
        
if __name__ == '__main__':
    
    exo = ExoCross('NaH')
    exo.check_input()
    # exo.load_PT_grid()
    exo.generate_inp(temperature=1000, pressure=1e-3)
    
    