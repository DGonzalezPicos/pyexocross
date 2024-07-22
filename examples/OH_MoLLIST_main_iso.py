from pyexocross.exocross import ExoCross

exo = ExoCross('OH_MoLLIST')

path_to_excross = '/data2/dario/pyexocross/exocross'
exo.set_path_exocross(path_to_excross)

exo.set_output('OH_MoLLIST_main_iso') # name of linelist
exo.load_PT_grid(file='PT_grids/PTpaths_high_4000.ls')
exo.xcross_grid(Nprocs=12, nice=False)

path_to_pRT = '/net/lem/data2/pRT_input_data/'
exo.set_path_pRT_input(path_to_pRT)
exo.copy_linelist_to_pRT()

# activate python environment
# $ source /data2/dario/retfish_env/activate.sh
# nohup python examples/OH_MoLLIST_main_iso.py > & OH_MoLLIST_run1.log &

# lastly, enable permissions to use linelists in LEM
# $ chmod -R 777 /net/lem/data2/pRT_input_data/OH_MoLLIST_main_iso