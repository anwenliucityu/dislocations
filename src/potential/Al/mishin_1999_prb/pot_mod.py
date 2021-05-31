#!/usr/bin/python3.8
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../..")))
from default_path import pot_path 

#########################################
#       information of potential        #
#    DP potential for Ti by Wenqi.T     #
#########################################

# lammps in file data
pot_name       =    'Ti_wenqi_2021'
pot_element    =    'Ti'
pot_type       =    'deepmd'
in_pot         =    ['pair_style         deepmd frozen_model.pb',
                     'pair_coeff', 
                     'neighbor           2.0 bin']
pot_path       =     [os.path.join(pot_path, 'potential/Ti/Ti_wenqi_2021/frozen_model.pb')]

# potential information
pot_cutoff     =    8.4

# information for slurm
module_load    =    'module load lammps/20201029/aocc_openmpi_deepmd_api_znver2'
appexe         =    'lmp'

