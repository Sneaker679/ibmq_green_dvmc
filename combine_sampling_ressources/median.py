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

# Fetching reference matrices
ref = []
ref.append(np.load(os.path.join(directories[0],os.path.join('output','H_AC_fock.npy'))))
ref.append(np.load(os.path.join(directories[0],os.path.join('output','H_CA_fock.npy'))))
ref.append(np.load(os.path.join(directories[0],os.path.join('output','S_AC_fock.npy'))))
ref.append(np.load(os.path.join(directories[0],os.path.join('output','S_CA_fock.npy'))))


final_shape = np.load(os.path.join(directories[0],os.path.join('output','H_AC.npy'))).shape
N = int(directory[0][0])

# Loading matrices
matrix_types = ['H_AC','H_CA','S_AC','S_CA']
matrices = [[],[],[],[]]
for index,matrix_type in enumerate(matrix_types):
    for index2, directory in enumerate(directories):
        directory_list_form = directory.split('_')
        directory_date,label = int(directory_list_form[-4]),int(directory_list_form[-2])

        matrix_path = os.path.join(directory,os.path.join('output',f'{matrix_type}.npy'))
        matrix = np.load(matrix_path)

        shape = matrix.shape

        matrices[index].append(matrix)
		

# Create 3D matrix
nd_H_AC = np.zeros((len(directories),final_shape[0],final_shape[1]))
nd_H_CA = np.zeros((len(directories),final_shape[0],final_shape[1]))
nd_S_AC = np.zeros((len(directories),final_shape[0],final_shape[1]))
nd_S_CA = np.zeros((len(directories),final_shape[0],final_shape[1]))
nd_matrices = [nd_H_AC,nd_H_CA,nd_S_AC,nd_S_CA]

for nd_index,nd_matrix in enumerate(nd_matrices):
    for index,matrix in enumerate(matrices[nd_index]):
        for row in range(final_shape[0]):
            for column in range(final_shape[1]):
                nd_matrices[nd_index][index,row,column] = matrix[row,column]


# Median
median_matrices = []
for nd_matrix in nd_matrices:
    median_matrices.append(np.median(nd_matrix,axis=0))

mu = int(directory_list_form[1].split('u')[1])
total_shots = sum(shots_list)
today = date.today()
date = today.strftime("%b %d %Y").split()
month,day = date[0].lower(),int(date[1])
label = 1
name_of_directory = f'median_{N}sites_mu{mu}_{day}_{month}_{label}_{total_shots}sh'

while os.path.exists(name_of_directory):
    label += 1
    name_of_directory = f'median_{N}sites_mu{mu}_{day}_{month}_{label}_{total_shots}sh'
    
os.mkdir(name_of_directory)
os.mkdir(os.path.join(name_of_directory,'output'))
for matrix in zip(median_matrices,matrix_types):
    np.save(os.path.join(name_of_directory,os.path.join('output',matrix[1])),matrix[0])

# Untested for not linux systems
os.system(f'cp files_to_copy/parameters.py {name_of_directory}')
os.system(f'cp files_to_copy/run_except_generate.sh {name_of_directory}')
