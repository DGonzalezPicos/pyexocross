from pyexocross.exocross import ExoCross
exo = ExoCross('NaH', path ='/data2/dario/pyexocross/')

path_to_excross = '/data2/dario/pyexocross/exocross'
exo.set_exocross_path(path_to_excross)

exo.check_input()
exo.read_definitions_file()
exo.set_output('NaH_rivlin') # name of linelist
exo.load_PT_grid(file='PT_grids/PTpaths_test.ls')
exo.xcross_grid(Nprocs=1)