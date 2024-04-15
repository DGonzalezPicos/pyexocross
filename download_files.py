'''Download .bz2 files from ExoMol using 'wget'.'''
import pathlib
# import glob
import subprocess
# import numpy as np

# Define the species and the url
# url = 'https://www.exomol.com/db/H2O/1H2-18O/HotWat78/'
# species = 'H2O_18'
# url = 'https://www.exomol.com/db/CN/12C-14N/Trihybrid/'
# url = 'https://www.exomol.com/db/CN/13C-14N/Trihybrid/'
url = 'https://www.exomol.com/db/H2O/1H2-18O/HotWat78/'

species = 'H2O_181'

# create folder
folder = pathlib.Path(f'{species}/input')
# folder = pathlib.Path(f'./')
folder.mkdir(parents=True, exist_ok=True)


# check for file containing the list of files to download
file = pathlib.Path(f'{species}/input/download_files.dat')
assert file.exists(), f'File {file} does not exist.'

with open(file, 'r') as f:
    lines = f.readlines()
    lines = [line.strip() for line in lines]

lines_bz2 = [line.split('.bz2')[0] + '.bz2' for line in lines if '.bz2' in line]

# save a log with the files to download and the files already downloaded
files_to_download =[f'{url}{line}' for line in lines_bz2]
# check for partition function files
lines_pf = [line.split('.pf')[0] + '.pf'for line in lines if '.pf' in line]
if len(lines_pf) > 0:
    files_to_download += [f'{url}{line}' for line in lines_pf]
    metadata_file = files_to_download[-1].replace('.pf', '.def')
    files_to_download += [metadata_file]

files_downloaded = sorted(folder.glob('*'))
files_downloaded = [file_i.name for file_i in files_downloaded]
# files_downloaded - []
print(f'Files to download: {files_to_download}')

print(f'\nFiles already downloaded: {files_downloaded}')

# check if partition function files are already downloaded
files_downloaded_pf = sorted(folder.glob('*.pf'))
assert len(files_downloaded_pf) <= 1, 'More than one partition function file found.'
if len(files_downloaded_pf) > 0:
    files_downloaded += files_downloaded_pf[0].name
print(f'Found {len(files_to_download)} files to download.')
print(f'Found {len(files_downloaded)} files already downloaded.')
download = True

if download:


    for i, file_i in enumerate(files_to_download):
        
        filename = file_i.split('/')[-1]
        # # remove extension if 
        filename = filename.split('.bz2')[0]
        # # add it back 
        # filename += '.bz2'
        
        if filename in files_downloaded:
            print(f'--> {filename} already downloaded.')
            # skip
            continue
       
        print(f'Downloading {i+1}/{len(file_i)}...', end='\r')
        subprocess.run(['wget', file_i], cwd=folder)

    # unzip all .bz2 files
    bz2_files = sorted(folder.glob('*.bz2'))
    for file_i in bz2_files:
        subprocess.run(['bzip2', '-d', file_i.name], cwd=folder)

    # copy this python file to the folder
    # subprocess.run(['cp', __file__, folder])
    # print(f'Copied {__file__} to {folder}')
    print('Done!')