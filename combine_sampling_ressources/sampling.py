import os,sys
from datetime import date
import numpy as np

np.set_printoptions(linewidth= 10000,precision=2,suppress=True)

# Fetching available directories
all_directories = [name for name in os.listdir(".") if os.path.isdir(name)]
available_directories = [x for x in all_directories
    if (not x[0] == 'c' 
    and not x[0] == 'f'
    and not x[0] == 't'
    and not x[0] == 'm'
    and not x[0] == 'e')]
available_directories.sort()

# Printing directories
print('Available directories:')
for index,dir in enumerate(available_directories):
    print(f'Index {index}  :  {dir}')
print('Input "done" to end the selection.')
directories = []
loop = True
while loop:
    answer = input('Directory to combine (list index): ')
    if answer.lower() == 'done':
        loop = False
    elif len(answer.split(':')) == 2:
        answer = answer.split(':')
        directories.extend(available_directories[int(answer[0]):int(answer[1])])
    else:
        answer = int(answer)
        directories.append(available_directories[answer])

# Determine shots
shots_list = []
for directory in directories:
    shot = directory.split('_')[-1][0:-2]
    shots_list.append(int(shot))

# Creating empty matrices of correct sizes.
matrices = []
matrix_types = ['H_AC','H_CA','S_AC','S_CA']
H_AC_sum = np.zeros(np.load(os.path.join(directories[0],os.path.join('output','H_AC.npy'))).shape)
H_CA_sum = np.zeros(np.load(os.path.join(directories[0],os.path.join('output','H_CA.npy'))).shape)
S_AC_sum = np.zeros(np.load(os.path.join(directories[0],os.path.join('output','S_AC.npy'))).shape)
S_CA_sum = np.zeros(np.load(os.path.join(directories[0],os.path.join('output','S_CA.npy'))).shape)
matrices.append(H_AC_sum)
matrices.append(H_CA_sum)
matrices.append(S_AC_sum)
matrices.append(S_CA_sum)
ref = []
ref.append(np.load(os.path.join(directories[0],os.path.join('output','H_AC_fock.npy'))))
ref.append(np.load(os.path.join(directories[0],os.path.join('output','H_CA_fock.npy'))))
ref.append(np.load(os.path.join(directories[0],os.path.join('output','S_AC_fock.npy'))))
ref.append(np.load(os.path.join(directories[0],os.path.join('output','S_CA_fock.npy'))))


print('\nPrint matrices?')
print('    Option 1: y')
print('    Option 2: n')
prt = input('Selected option: ')
if prt.lower() == 'y':
    prt = True
else:
    prt = False

# Appliying modifications to the matrices, and combining them using an average.
N = int(directories[0][0])
for index,matrix_type in enumerate(matrix_types):
    for ii in zip(directories,shots_list):
        directory = ii[0]
        directory_list_form = directory.split('_')
        directory_date,label = int(directory_list_form[-4]),int(directory_list_form[-2])

        matrix_path = os.path.join(ii[0],os.path.join('output',f'{matrix_type}.npy'))
        matrix = np.load(matrix_path)

        matrices[index] = np.add(matrices[index],matrix*(ii[1])/sum(shots_list))

        shape = matrix.shape

        if prt is True:
            print(matrix_type,directory)
            zeros = np.zeros((shape[0],1))
            print(np.concatenate((matrix,zeros,matrices[index],zeros,ref[index]),axis=1))
            print()

#name_of_directory = input('Name of the new directory: ')
mu = int(directory_list_form[1].split('u')[1])
total_shots = sum(shots_list)
today = date.today()
date = today.strftime("%b %d %Y").split()
month,day = date[0].lower(),int(date[1])
label = 1
name_of_directory = f'combined_{N}sites_mu{mu}_{day}_{month}_{label}_{total_shots}sh'

while os.path.exists(name_of_directory):
    label += 1
    name_of_directory = f'combined_{N}sites_mu{mu}_{day}_{month}_{label}_{total_shots}sh'
    
os.mkdir(name_of_directory)
os.mkdir(os.path.join(name_of_directory,'output'))
for matrix in zip(matrices,matrix_types):
    np.save(os.path.join(name_of_directory,os.path.join('output',matrix[1])),matrix[0])

# Untested for not linux systems
os.system(f'cp files_to_copy/parameters.py {name_of_directory}')
os.system(f'cp files_to_copy/run_except_generate.sh {name_of_directory}')
