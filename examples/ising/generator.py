# Этот код генерирует систему Изинга размером L*L, создает ini-файл для запуска через метрополис 
# и выдает в консоль команду для запуска программы на кластере 
# За основу берется файл sample.ini
# Запустить в цикле эту программу можно так:
# for i in {10..32..2}; do python generator.py $i; done;


import numpy as np
import matplotlib.pyplot as plt
import sys
import configparser

if len(sys.argv) != 2:
    raise AttributeError("run format: %s <L>"%(sys.argv[0]))

L = int(sys.argv[1])
target_mfsys = "ising_%d.mfsys"%(L*L)
target_ini = "ising_%d.ini"%(L*L)
source_ini = "sample.ini"

def chess_board(n):
    # Индексы строк и столбцов
    rows = np.arange(n)[:, np.newaxis]
    cols = np.arange(n)
    
    # (i + j) % 2 определяет цвет клетки
    board = np.where((rows + cols) % 2 == 0, 1, -1)
    return board

def save_to_mfsys(filename, system):
    N = system.shape[1]
    head = f'[header]\nversion=2\ndimensions=3\ntype=standart\nsize={N}\nemin=0.000000\nemax=0.000000\nstate={"0"*N}\nminstate=\nmaxstate=\ninteractionrange=4.5\nsizescale=1\nmagnetizationscale=1\n[parts]'
    np.savetxt(filename,
               np.vstack((
                   np.arange(N), 
                   system[0], 
                   system[1],
                   np.zeros(N), 
                   np.zeros(N), 
                   np.zeros(N), 
                   system[2], 
                   np.zeros(N))).T, 
                fmt="%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i", comments="", header=head)


X,Y = np.meshgrid(np.arange(L),np.arange(L))
Mz = chess_board(L)
sys = np.vstack((X.reshape(-1),Y.reshape(-1),Mz.reshape(-1)))
save_to_mfsys(target_mfsys,sys)