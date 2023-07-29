import os,sys
import numpy as np

np.set_printoptions(linewidth= 10000,precision=2,suppress=True)

available_directories = [name for name in os.listdir(".") if os.path.isdir(name)]
available_directories = [x for x in available_directories
    if (not x[0] == 'c' 
    and not x[0] == 'f'
    and not x[0] == 'e')]
available_directories.sort()
print('Available directories:')
for index,dir in enumerate(available_directories):
    print(f'Index {index}  :  {dir}')
print('Input "done" to end the selection.')
directories = []
finish = True
while finish:
    answer = input('Directory to combine (list index): ')
    if answer.lower() == 'done':
        finish = False
    else:
        answer = int(answer)
        directories.append(available_directories[answer])

# Determine shots
shots_list = []
for directory in directories:
    shot = directory.split('_')[-1][0:-2]
    shots_list.append(int(shot))

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


# Load matrices
for index,matrix_type in enumerate(matrix_types):
    for ii in zip(directories,shots_list):
        directory = ii[0]
        directory_date,label = int(directory.split('_')[-4]),int(directory.split('_')[-2])

        matrix_path = os.path.join(ii[0],os.path.join('output',f'{matrix_type}.npy'))
        matrix = np.load(matrix_path)

        matrices[index] = np.add(matrices[index],matrix*(ii[1])/sum(shots_list))

name_of_directory = input('Name of the new directory: ')
os.mkdir(name_of_directory)
os.mkdir(os.path.join(name_of_directory,'output'))
for matrix in zip(matrices,matrix_types):
    np.save(os.path.join(name_of_directory,os.path.join('output',matrix[1])),matrix[0])

os.system(f'cp files_to_copy/parameters.py {name_of_directory}')
os.system(f'cp files_to_copy/run_except_generate.sh {name_of_directory}')
