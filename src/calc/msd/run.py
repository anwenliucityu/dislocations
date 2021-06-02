#!/usr/bin/python3.8
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),".."))) #yes
import default_path as path
spec      = __import__(sys.argv[2].replace('.py',''))
pot_path  = path.pot_path
sys.path.append(os.path.abspath(os.path.join(pot_path,'potential',spec.pot_element,spec.pot_id)))
import pot_mod
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../dislocation_boy/msd")))
import msd
exec(f'import info_{spec.element_struct}_{spec.start_temperature}K as mi')
import importlib

# Read information in spec file and path
main_path = path.output_path

#####################
#   Main Function   #
#####################

# Customize
params = {}

if 'sbatch_job' in dir(spec):
    params['sbatch_job']            = spec.sbatch_job
if 'ovito' in dir(spec):
    params['ovito']                 = spec.ovito

unit_cell_size = [spec.num_unit_cell_x, spec.num_unit_cell_y, spec.num_unit_cell_z]

if spec.element_struct =='fcc':
    config = 'fcc_screw'
if spec.element_struct =='hcp':
    config = 'hcp_screw_a_basal'
if spec.element_struct =='bcc':
    config = 'bcc_screw'

# import the module to constructure the final atomic structure
module_path = os.path.abspath(msd.__file__)

# change module import and reimport
lines = open(module_path, 'r').readlines()
line_zero = lines.pop(3)
line_zero = line_zero.replace(line_zero.split()[3], config)
lines.insert(3, line_zero)
line_zero = lines.pop(7)
line_zero = line_zero.replace(line_zero.split()[3], config)
lines.insert(7, line_zero)
open(module_path,'w').writelines(lines)
importlib.reload(msd)

if __name__ == "__main__":
    
    if sys.argv[1] == 'calc':
        msd.calc_msd(main_path, spec.pot_element, spec.pot_id, unit_cell_size, spec.element_struct,
                    spec.start_temperature, mi.latt_const, mi.mass, pot_mod.pot_path, 
                    spec.warm_runningstep, spec.ave_interval, spec.ave_times, pot_mod.in_pot,
                    spec.partition, spec.mem, pot_mod.module_load, pot_mod.appexe, spec.ncpu,
                    **params)

    if sys.argv[2] ==  'plot':

        a =  True
