import os
import numpy as np
from voigt_notation import lmp_stress_atom

def output_dir(main_path, pot_element, pot_name, element_struct, config_style, dislocation_type, 
               simulation_type, unit_cell_size ,global_emin=True, T=0, temp=300):
    temp = str(int(temp)) + 'K'
    T    = str(int(T)) + 'K'
    box  = str(unit_cell_size[0])+'_'+str(unit_cell_size[1])+'_'+ str(unit_cell_size[2])
    dir = os.path.join(main_path, pot_element, pot_name, simulation_type, element_struct, config_style, dislocation_type)
                 # e.g. main/Ti/Ti_kevin_2020/energy_minimization/hcp/cylinder/screw_a_basal/
    if config_style == 'quadrupolar':
        global_emin = False
    if simulation_type == 'energy_minimization':
        if global_emin == True:
            dir = os.path.join(dir, T+'_'+temp+'_'+T,)
            # e.g. energy_minimization/0_300_0
        if global_emin == False:
            dir = os.path.join(dir, T)
    elif simulation_type == 'metastable':
        dir = os.path.join(dir, T+'_'+temp)
    elif simulation_type == 'finite_T':
        dir = os.path.join(dir, T)
    dir = os.path.join(dir, box)
    if not os.path.isdir(dir):
        os.makedirs(dir)
    dump_dir = os.path.join(dir, 'dump_file')
    if not os.path.isdir(dump_dir):
        os.makedirs(dump_dir)
    print("Atomic structure is prepared ......")
    print("work_path = " + dir)
    return dir

def write_cfg_auxiliary(atom_coor, pot_element, dislocation_type, mass, box_boundary, directory,
                     **kwargs):
    '''
    write cfg with auxiliary properties
    '''

    # size of box (Angstrom)
    box_size = [box_boundary[0][1] - box_boundary[0][0],
                box_boundary[1][1] - box_boundary[1][0],
                box_boundary[2][1] - box_boundary[2][0]]

    # transform fractional coordinates
    num_atom = atom_coor.shape[0]
    atom_fractional_coor = np.empty(shape=(num_atom, 3))
    for axis in range(3):
        atom_fractional_coor[:,axis] = \
                (atom_coor[:,axis] - box_boundary[axis][0]) / box_size[axis]

    case_name = pot_element + '_' + dislocation_type
    filename = case_name + '.cfg'
    file_path_name = os.path.join(directory, filename)

    file = open(file_path_name, 'w')

    file.write(f"Number of particles = {num_atom}\n"
                "A = 1.0 Angstrom (basic length-scale)\n"
               f"H0(1,1) = {box_size[0]} A\n"
                "H0(1,2) = 0 A\n"
                "H0(1,3) = 0 A\n"
                "H0(2,1) = 0 A\n"
               f"H0(2,2) = {box_size[1]} A\n"
                "H0(2,3) = 0 A\n"
                "H0(3,1) = 0 A\n"
                "H0(3,2) = 0 A\n"
               f"H0(3,3) = {box_size[2]} A\n"
                ".NO_VELOCITY.\n"
               f"entry_count = {3 + len(kwargs)}\n")
    i = 0
    for key, _ in kwargs.items():
        file.write(f"auxiliary[{i}] = {key}\n")
        i += 1

    file.write(f"{mass}\n"
              f"{pot_element}\n")

    for i in range(num_atom):
        file.write(f"{atom_fractional_coor[i,0]:>14.15f} "
                   f"{atom_fractional_coor[i,1]:>14.15f} "
                   f"{atom_fractional_coor[i,2]:>14.15f} ")
        for _, value in kwargs.items():
            file.write(f"{value[i]} ")

        file.write("\n")

    file.close()
 
    
def write_datafile(pot_element, dislocation_type, mass, directory, atom_coor, atom_type, box_boundary, suffix='_disl.dat'):
    '''write datafile'''

    num_atom = atom_coor.shape[0]
    type_list = list(np.unique(atom_type))
    type_list.sort()
    num_type = len(type_list)

    case_name = pot_element + '_' + dislocation_type
    filename = case_name + suffix
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"{case_name}\n\n"
               f"{num_atom} atoms\n"
               f"{num_type} atom types\n"
               f"{box_boundary[0][0]:>28.20f} {box_boundary[0][1]:>28.20f} xlo xhi\n"
               f"{box_boundary[1][0]:>28.20f} {box_boundary[1][1]:>28.20f} ylo yhi\n"
               f"{box_boundary[2][0]:>28.20f} {box_boundary[2][1]:>28.20f} zlo zhi\n\n")

    
    file.write("Masses\n\n")
    for each_type in type_list:
        file.write(f"{int(each_type):2d} {mass}\n")
    file.write("\n")

    file.write("Atoms\n\n")
    for i in range(num_atom):
        file.write(f"{i+1:>8d} {int(atom_type[i]):>4d} \
                   {atom_coor[i,0]:>28.20f} {atom_coor[i,1]:>28.20f} {atom_coor[i,2]:>28.20f}\n")

    file.close()
    
def write_in_file(pot_path, config_style, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type, 
                  mass, shell_radius, sample_center, disl_center, calc_atomic_stress=False, global_emin=False, temp = 300.0,
                  cooling_rate=500):
  
    if element_struct == 'fcc' :
        atom_vol = str(latt_const**3/4/100000)
    if element_struct == 'hcp' :
        atom_vol = str(np.sqrt(3)*latt_const["a"]**2*latt_const["c"]/4/100000)
    elem_symbol = pot_element
    case_name = pot_element + '_' + dislocation_type
    calc = "_* mass type x y z\n"

    ini_path = os.getcwd()
    os.chdir(directory)
    for i in range(len(pot_path)):
        pot_name = os.path.basename(pot_path[i])
        os.system(f'ln -s {pot_path[i]} {pot_name}')
    os.chdir(ini_path)
 
    filename = 'in.lammps'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')
    file.write("#----------------------------------------#\n"
               "#---   basic info & initialization    ---#\n"
               "#----------------------------------------#\n")
    file.write(f"clear\n"
                "units         metal\n"
                "dimension     3\n"
                "boundary      f f p\n"
                "atom_style    atomic\n"
                "atom_modify   map array\n"
               f"read_data     {case_name}_disl.dat\n"
               f"mass          1 {mass}\n\n") 
    
    file.write("#----------------------------------------#\n"
               "#---            potential             ---#\n"
               "#----------------------------------------#\n")

    for i in range(len(in_pot)):
        file.write(f"{in_pot[i]}\n")
    file.write('\n')

    if config_style != 'quadrupolar':
        file.write("#----------------------------------------#\n"
                   "#---     freeze atoms in boundary     ---#\n"
                   "#----------------------------------------#\n")

        file.write(f"region    inner cylinder z {sample_center[0]} {sample_center[1]} {shell_radius} EDGE EDGE\n"
                    "group     free_atom region inner\n"
                    "group     freeze_atom subtract all free_atom\n"
                    "fix       freeze_boundary freeze_atom setforce 0.0 0.0 0.0\n\n")
    file.write("shell     cd dump_file\n\n")

    if global_emin == True:
        running_step = temp * cooling_rate 
        file.write("#----------------------------------------#\n"
                   "#---         heat to help Emin        ---#\n"
                   "#----------------------------------------#\n")
        file.write("compute        freeatomtemp all temp\n"
                   "thermo         20\n"
                   "thermo_style   custom step pe lx ly lz press c_freeatomtemp\n"
                  f"velocity       free_atom create {temp} 4928459 rot yes dist gaussian\n"
                  f"fix            warm free_atom nvt temp {temp} {temp} 0.1\n"
                  f"run            2000\n"
                   "unfix          warm\n"
                  f"fix            cooling free_atom nvt temp {temp} 0.1 0.1\n"
                  f"run            {running_step}\n"
                   "unfix          cooling\n\n")
    
    if calc_atomic_stress == True:
        calc = "_* mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
        file.write("#----------------------------------------#\n"
                   "#---              compute             ---#\n"
                   "#----------------------------------------#\n")
        file.write("compute   peratom all stress/atom NULL\n")
        for i in range(1,7):
            file.write(f"variable    s{lmp_stress_atom[i]} atom c_peratom[{i}]/" + atom_vol + '\n')
        file.write("\n")

    file.write("#----------------------------------------#\n"
               "#---                dump              ---#\n"
               "#----------------------------------------#\n")
  
    file.write(f"dump            emin all custom 200 dump.{case_name}{calc}" 
                "dump_modify     emin element Ti sort id\n\n"
                "reset_timestep  0\n"
                "thermo          50\n"
                "thermo_style    custom step pe lx ly lz press\n\n")

    file.write("#----------------------------------------#\n"
               "#---         energy minimization      ---#\n"
               "#----------------------------------------#\n"
               "min_style cg\n"
               "minimize 1e-25 1e-25 5000 10000\n\n"
               "print"' "'"All done"'"' )
  
    file.close()


def write_in_finite_T_meta_file(pot_path, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, shell_radius, sample_center, disl_center, simulation_type, start_temp,
                  spring_factor_k = None,
                  running_steps=None, dump_interval=None, calc_atomic_stress=False, temp = 300.0):

    if element_struct == 'fcc' :
        atom_vol = str(latt_const**3/4/100000)
    if element_struct == 'hcp' :
        atom_vol = str(np.sqrt(3)*latt_const["a"]**2*latt_const["c"]/4/100000)
    elem_symbol = pot_element
    case_name = pot_element + '_' + dislocation_type
    calc = "_* mass type x y z\n"

    ini_path = os.getcwd()
    os.chdir(directory)
    for i in range(len(pot_path)):
        pot_name = os.path.basename(pot_path[i])
        os.system(f'ln -s {pot_path[i]} {pot_name}')
    os.chdir(ini_path)

    filename = 'in.lammps'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')
    file.write("#----------------------------------------#\n"
               "#---   basic info & initialization    ---#\n"
               "#----------------------------------------#\n")
    file.write(f"clear\n"
                "units         metal\n"
                "dimension     3\n"
                "boundary      f f p\n"
                "atom_style    atomic\n"
                "atom_modify   map array\n"
               f"read_data     {case_name}_disl.dat\n"
               f"mass          1 {mass}\n\n")
    
    file.write("#----------------------------------------#\n"
               "#---            potential             ---#\n"
               "#----------------------------------------#\n")

    for i in range(len(in_pot)):
        file.write(f"{in_pot[i]}\n")

    file.write("\n"
               "#----------------------------------------#\n"
               "#---     freeze atoms in boundary     ---#\n"
               "#----------------------------------------#\n")

    file.write(f"region    inner cylinder z {sample_center[0]} {sample_center[1]} {shell_radius} EDGE EDGE\n"
                "group     free_atom region inner\n"
                "group     freeze_atom subtract all free_atom\n")
    if simulation_type == 'metastable':
        file.write("fix       freeze_boundary freeze_atom setforce 0.0 0.0 0.0\n"
                   "shell     cd dump_file\n\n")
    if simulation_type == 'finite_T':
        temp = start_temp
        file.write(f"fix       freeze_boundary freeze_atom spring/self {spring_factor_k}\n"
                    "shell     cd dump_file\n\n")

    file.write("#----------------------------------------#\n"
               "#---          add temperature         ---#\n"
               "#----------------------------------------#\n")                                                                                
    file.write("compute        freeatomtemp all temp\n"
               "thermo         100\n"
               "thermo_style   custom step pe lx ly lz press c_freeatomtemp\n"
              f"velocity       free_atom create {temp} 4928459 rot yes dist gaussian\n"
              f"fix            warm free_atom nvt temp {temp} {temp} 0.1\n"
              f"run            3000\n"
               "reset_timesetp 0\n\n")

    if calc_atomic_stress == True:
        calc = "_* mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
        file.write("#----------------------------------------#\n"
                   "#---              compute             ---#\n"
                   "#----------------------------------------#\n")
        file.write("compute   peratom all stress/atom NULL\n")
        for i in range(1,7):
            file.write(f"variable    s{lmp_stress_atom[i]} atom c_peratom[{i}]/" + atom_vol + '\n')
        file.write("\n")

    file.write("#----------------------------------------#\n"
               "#---                dump              ---#\n"
               "#----------------------------------------#\n")
  
    file.write(f"dump           {simulation_type} all custom {dump_interval} dump.{case_name}{calc}" 
               f"dump_modify    {simulation_type} element Ti sort id\n\n"
                "reset_timestep 0\n"
                "thermo         200\n"
                "thermo_style   custom step pe lx ly lz press\n"
               f"run            {running_steps}\n\n")

    file.write("print"' "'"All done"'"' )
  
    file.close()



  
def write_job(pot_element, dislocation_type, partition, mem, module_load, appexe, directory, ncpu=12):
    '''
    write lammps in file
    '''
    case_name = pot_element + '_' + dislocation_type
    in_filename = 'in.lammps'
    core_num = str(ncpu) 
   
    filename = 'job.sh'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')
  
    file.write(f"#!/usr/bin/env bash\n\n"
                "#SBATCH --job-name="'"' + case_name + '"\n'
               f"#SBATCH --partition={partition}\n"
                "#SBATCH --ntasks-per-core=1\n"
               f"#SBATCH --ntasks={core_num}\n"
                "#SBATCH --nodes=1\n"
                "#SBATCH --cpus-per-task=1\n"
               f"#SBATCH --mem={mem}\n"
                "#SBATCH --error='"'j%j.stderr'"'\n"
                "#SBATCH --output='"'j%j.stdout'"'\n\n"
                "# import some basic util func\n"
                ". $HOME/template/shell/bash/bash_script_home\n"
                ". $HOME/template/shell/bash/sys_util\n\n"
                "export OMP_NUM_THREADS=$'{'SLURM_CPUS_PER_TASK'}'\n"
                )

    file.write(f"ulimit -s unlimited\n\n"
               f"{module_load}\n"
                "appexe='"f'nice -n20 {appexe}'"'\n\n"
               f"input_file={in_filename}\n"
                "output_file=$input_file.stdout\n"
                "log_file=$input_file.log\n\n"
                "mpirun -np $SLURM_NTASKS $appexe -in $input_file -log $log_file > $output_file\n")
    file.close()
    
    
