#!/usr/bin/python3.8
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../.."))) #yes
import default_path as path
spec      = __import__(sys.argv[2].replace('.py',''))
pot_path  = path.pot_path
sys.path.append(os.path.abspath(os.path.join(pot_path,'potential',spec.pot_element,spec.pot_id)))
import pot_mod
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),"../../../dislocation_boy/disl_aniso")))
#import material_info as mi
exec(f'import info_{spec.element_struct}_{spec.start_temperature}K as mi')
import importlib

# Read information in spec file and path
main_path = path.output_path

#####################
#   Main Function   #
#####################

# Customize
params = {}
if 'boundary_freeze_width' in dir(spec):
    params['boundary_freeze_width'] = spec.boundary_freeze_width
if 'calc_atomic_stress' in dir(spec):
    params['calc_atomic_stress']    = spec.calc_atomic_stress
if 'global_emin' in dir(spec):
    params['global_emin']           = spec.global_emin
if 'cooling_rate' in dir(spec):
    params['cooling_rate']          = spec.cooling_rate
if 'calc_aniso_stress' in dir(spec):
    params['calc_aniso_stress']     = spec.calc_aniso_stress
if 'temp' in dir(spec):
    params['temp']                  = spec.temp
if 'running_steps' in dir(spec):
    params['running_steps']         = spec.running_steps
if 'dump_interval' in dir(spec):
    params['dump_interval']         = spec.dump_interval
if 'sbatch_job' in dir(spec):
    params['sbatch_job']            = spec.sbatch_job
if 'spring_factor_k' in dir(spec):
    params['spring_factor_k']       = spec.spring_factor_k
if 'output_perfect' in dir(spec):
    params['output_perfect']        = spec.output_perfect
if 'ovito' in dir(spec):
    params['ovito']                 = spec.ovito

unit_cell_size = [spec.num_unit_cell_x, spec.num_unit_cell_y, spec.num_unit_cell_z]
disl_center    = [spec.disl_center_x, spec.disl_center_y]

# import the module to constructure the final atomic structure
if spec.config_style == 'quadrupolar':
    import quadru_aniso_disl_config_constructor as aniso_disl
else:
    import aniso_disl_config_constructor as aniso_disl
module_path = os.path.abspath(aniso_disl.__file__)

# change module import and reimport
lines = open(module_path, 'r').readlines()
line_zero = lines.pop(0)
line_zero = line_zero.replace(line_zero.split()[3], spec.element_struct+'_'+spec.dislocation_type)
open(module_path,'w').writelines(line_zero)
open(module_path,'a').writelines(lines)
importlib.reload(aniso_disl)

if __name__ == "__main__":
    
    if sys.argv[1] == 'calc':
        result = aniso_disl.aniso_disl_constructor(main_path, spec.config_style, spec.dislocation_type, 
                           spec.pot_element, mi.latt_const, spec.pot_id, spec.element_struct, unit_cell_size,
                           spec.simulation_type, spec.start_temperature,
                           mi.elastic_const, mi.mass, pot_mod.pot_cutoff, pot_mod.in_pot, pot_mod.pot_type, pot_mod.pot_path,
                           disl_center,
                           spec.partition, spec.mem, pot_mod.module_load, pot_mod.appexe, spec.ncpu,
                           **params)
    if sys.argv[2] ==  'plot':

        a =  True
