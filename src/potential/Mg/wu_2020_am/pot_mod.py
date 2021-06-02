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
pot_element    =    'Mg'
pot_type       =    'meam/c'
in_pot         =    ['pair_style         meam/c',
                     'pair_coeff         * * library.meam Mg Mg.meam Mg',
                     'neighbor           0.3 bin']
pot_path       =     [os.path.join(pot_path, 'potential/Mg/wu_2021_am/library.meam'),
                      os.path.join(pot_path, 'potential/Mg/wu_2021_am/Mg.meam')]

# potential information
pot_cutoff     =    8.4

# information for slurm
module_load    =    'module load lammps/gcc_openmpi_zen2'
appexe         =    'lmp'

