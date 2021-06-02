import os
import numpy as np
from elem_dict import elem_dict
from voigt_notation import lmp_stress_atom

def output_directory(config, elem_property, directory_path, ):
    '''
    directory_path has the format of a directory path
    '''
    element = elem_property["elem_symbol"]
    case_name = config["case_name"]
    cell_x_number = '_size_' + str(config["box_num_cell_x"]) + '_'
    cell_y_number = str(config["box_num_cell_y"]) + '_'
    cell_z_number = str(config["box_num_cell_z"]) + '/'
    directory = directory_path
    directory += '/' + element
    directory += '/' + case_name + cell_x_number + cell_y_number + cell_z_number
    if not os.path.isdir(directory):
      os.makedirs(directory)  
    print("Please cd " + directory + "to see lammps calculated results")
    directory_dump = directory + '/dump_file'
    if not os.path.isdir(directory_dump):
      os.mkdir(directory_dump)
    return directory

def write_info(config, elem_property, atom_coor, atom_index_selected, sample_radius, shell_radius, box_boundary, directory):
    atom_num = str(atom_coor.shape[0])
    freeze_atom = str(len(atom_index_selected))
    free_atom = str(atom_coor.shape[0]-len(atom_index_selected))
    el_con = np.array(elem_property["elastic_const_initial"])/10e9
    case_name = config["case_name"]
    filename = 'basic_info.txt'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')
  
    file.write(f"Please NOTE that before sbatch job.sh file, please assign \n"
                "a gauss and core number for it to work, default is gauss12\n"
                "with 32 cores\n\n"
                "Configuration type           = " + config["case_name"] + "\n"
                "lattice constant a (A)       = " + str(elem_property["latt_const"]) + "\n"
               )
    
    file.write(f"elastic constant (GPa)       = |{el_con[0]:>7.1f}{el_con[1]:>7.1f}{el_con[2]:>7.1f}{el_con[3]:>7.1f}{el_con[4]:>7.1f}{el_con[5]:>7.1f}   |\n"
               f"                               |{el_con[6]:>14.1f}{el_con[7]:>7.1f}{el_con[8]:>7.1f}{el_con[9]:>7.1f}{el_con[10]:>7.1f}   |\n"
               f"                               |{el_con[11]:>21.1f}{el_con[12]:>7.1f}{el_con[13]:>7.1f}{el_con[14]:>7.1f}   |\n"
               f"                               |{el_con[15]:>28.1f}{el_con[16]:>7.1f}{el_con[17]:>7.1f}   |\n"
               f"                               |{el_con[18]:>35.1f}{el_con[19]:>7.1f}   |\n"
               f"                               |{el_con[20]:>42.1f}   |\n"
                )
                     
    file.write(f"Total simulated atom number  = " + atom_num + "\n"
                "Freeze atom number           = " + freeze_atom + "\n"
                "Free atom number             = " + free_atom + "\n\n"
                "Sample radius                = " + str(sample_radius) + " Angstrom\n"
                "Free-relax radius            = " + str(shell_radius) + " Angstrom\n\n"
                "number of repeat cell in X   = " + str(config["box_num_cell_x"]) + "\n"
                "number of repeat cell in Y   = " + str(config["box_num_cell_y"]) + "\n"
                "number of repeat cell in Z   = " + str(config["box_num_cell_z"]) + "\n"
                "box length in X              = " + str(-box_boundary[0][0]+box_boundary[0][1]) + " Angstrom\n"
                "box length in Y              = " + str(-box_boundary[1][0]+box_boundary[1][1]) + " Angstrom\n"
                "box length in Z              = " + str(-box_boundary[2][0]+box_boundary[2][1]) + " Angstrom\n\n"
                "box X direction              = " + str(np.transpose(config["frame_new"])[0]) + "\n"
                "box Y direction              = " + str(np.transpose(config["frame_new"])[1]) + "\n"
                "box Z direction              = " + str(np.transpose(config["frame_new"])[2]) + "\n"
  )
  
    file.close()

def write_cfg_auxiliary(atom_coor, config, elem_property, box_boundary, directory,
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
    
    case_name = config["case_name"]
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

    file.write(f"{elem_property['mass']}\n"
              f"{elem_property['elem_symbol']}\n")

    for i in range(num_atom):
        file.write(f"{atom_fractional_coor[i,0]:>14.15f} "
                   f"{atom_fractional_coor[i,1]:>14.15f} "
                   f"{atom_fractional_coor[i,2]:>14.15f} ")
        for _, value in kwargs.items():
            file.write(f"{value[i]} ")

        file.write("\n")

    file.close()
 
    
def write_datafile(config, directory, atom_coor, atom_type, box_boundary, tilt_para=[0.0,0.0,0.0]):
    '''write datafile'''

    if_write_mass = False

    num_atom = atom_coor.shape[0]
    type_list = list(set(atom_type))
    type_list.sort()
    num_type = len(type_list)

    case_name = config["case_name"]
    filename = case_name + '_aniso.dat'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"{case_name}\n\n"
               f"{num_atom} atoms\n"
               f"{num_type} atom types\n"
               f"{box_boundary[0][0]:>28.20f} {box_boundary[0][1]:>28.20f} xlo xhi\n"
               f"{box_boundary[1][0]:>28.20f} {box_boundary[1][1]:>28.20f} ylo yhi\n"
               f"{box_boundary[2][0]:>28.20f} {box_boundary[2][1]:>28.20f} zlo zhi\n"
               f"{tilt_para[0]:>28.20f} {tilt_para[1]:>28.20f} {tilt_para[2]} xy xz yz\n\n")

    if if_write_mass:
        file.write("Masses\n\n")
        for each_type in type_list:
            mass = elem_dict[config["type_to_elem"][each_type]]["mass"]
            file.write(f"{int(each_type):2d} {mass}\n")
        file.write("\n")

    file.write("Atoms\n\n")
    for i in range(num_atom):
        file.write(f"{i+1:>8d} {int(atom_type[i]):>4d} \
                   {atom_coor[i,0]:>28.20f} {atom_coor[i,1]:>28.20f} {atom_coor[i,2]:>28.20f}\n")

    file.close()
   
def write_perfect_datafile(config, directory, atom_coor, atom_type, box_boundary):
    '''write datafile'''

    if_write_mass = False

    num_atom = atom_coor.shape[0]
    type_list = list(set(atom_type))
    type_list.sort()
    num_type = len(type_list)

    case_name = config["case_name"]
    filename = case_name + '_perfect_reference.dat'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"{case_name}\n\n"
               f"{num_atom} atoms\n"
               f"{num_type} atom types\n"
               f"{box_boundary[0][0]:>28.20f} {box_boundary[0][1]:>28.20f} xlo xhi\n"
               f"{box_boundary[1][0]:>28.20f} {box_boundary[1][1]:>28.20f} ylo yhi\n"
               f"{box_boundary[2][0]:>28.20f} {box_boundary[2][1]:>28.20f} zlo zhi\n\n")

    if if_write_mass:
        file.write("Masses\n\n")
        for each_type in type_list:
            mass = elem_dict[config["type_to_elem"][each_type]]["mass"]
            file.write(f"{int(each_type):2d} {mass}\n")
        file.write("\n")

    file.write("Atoms\n\n")
    for i in range(num_atom):
        file.write(f"{i+1:>8d} {int(atom_type[i]):>4d} \
                   {atom_coor[i,0]:>28.20f} {atom_coor[i,1]:>28.20f} {atom_coor[i,2]:>28.20f}\n")

    file.close()
    
def write_in_file(config, elem_property, directory, shell_radius, sample_center, dislocation_line_position):
  
    pot_path = "/gauss12/home/cityu/anwenliu/git/anisotropic-elasticity/src/anwen_disl_aniso/potential/"
    if elem_property["cryst_structure"] == 'fcc' :
        atom_vol = str(elem_property["latt_const"]**3/4/100000)
    if elem_property["cryst_structure"] == 'hcp' :
        atom_vol = str(np.sqrt(3)*elem_property["latt_const"]["a"]**2*elem_property["latt_const"]["c"]/4/100000)
    elem_symbol = elem_property["elem_symbol"]
    case_name = config["case_name"]
    pot_file_name = elem_property["pot_file_name"]
    mass = str(elem_property["mass"])
 
    filename = 'in.'+ case_name
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"clear\n"
                "units metal\n"
                "dimension 3\n"
                "boundary f f p\n"
                "atom_style atomic\n"
                "atom_modify map array\n\n"
                "read_data " + case_name + "_aniso.dat\n\n"
                "mass 1 " + mass + "\n\n") 
    
    file.write(f"region inner cylinder z {sample_center[0]} {sample_center[1]} {shell_radius} EDGE EDGE\n"
                "group free_atom region inner\n"
                "group freeze_atom subtract all free_atom\n"
                "fix freeze_boundary freeze_atom setforce 0.0 0.0 0.0\n\n"
                "pair_style " + elem_property["pair_style"] + '\n'
                "pair_coeff * * ")
    for i in range(len(pot_file_name)):
        file.write(pot_path + pot_file_name[i] + ' ' + elem_symbol + ' ')
      
    file.write( "\n"
                "neighbor 2.0 bin\n\n"
                "compute peratom all stress/atom NULL\n\n")
  
    for i in range(1,7):
      file.write(f"variable s{lmp_stress_atom[i]} atom c_peratom[{i}]/" + atom_vol + '\n')
    file.write(f"\n")
  
    file.write(f"shell cd dump_file\n"
                "dump 11 all custom 500 dump." + case_name + "_*.cfg mass type xs ys zs v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n\n"
                "reset_timestep 0\n"
                "thermo 10\n"
                "thermo_style custom step pe lx ly lz press\n"
                "min_style cg\n"
                "minimize 1e-25 1e-25 5000 10000\n\n"
                "shell cd ../../\n"
                "shell mkdir result\n"
                "shell cd result\n"
                "dump dump_final all custom 1 "
                "dump_final_" + str(dislocation_line_position[0]*100) + str(dislocation_line_position[1]*100) +
                ".cfg mass type xs ys zs v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
                "run 0\n"
                "undump dump_final\n"
                "write_restart 0K_Emin.restart\n"
                "print"' "'"All done"'"' )
  
    file.close()

def write_dp_in_file(config, elem_property, directory, shell_radius, sample_center, dislocation_line_position):

    atom_vol = str(np.sqrt(3)*elem_property["latt_const"]["a"]**2*elem_property["latt_const"]["c"]/4/100000)
    elem_symbol = elem_property["elem_symbol"]
    case_name = config["case_name"]
    pot_file_name = elem_property["pot_file_name"]
    mass = str(elem_property["mass"])

    filename = 'in.'+ case_name
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"clear\n"
                "units metal\n"
                "dimension 3\n"
                "boundary f f p\n"
                "atom_style atomic\n"
                "atom_modify map array\n\n"
                "read_data " + case_name + "_aniso.dat\n\n"
                "mass 1 " + mass + "\n\n")

    file.write(f"region inner cylinder z {sample_center[0]} {sample_center[1]} {shell_radius} EDGE EDGE\n"
                "group free_atom region inner\n"
                "group freeze_atom subtract all free_atom\n"
                "fix freeze_boundary freeze_atom setforce 0.0 0.0 0.0\n\n"
                "pair_style deepmd frozen_model.pb \n"
                "pair_coeff    ")

    file.write( "\n"
                "neighbor 2.0 bin\n\n"
                "compute peratom all stress/atom NULL\n\n")

    for i in range(1,7):
      file.write(f"variable s{lmp_stress_atom[i]} atom c_peratom[{i}]/" + atom_vol + '\n')
    file.write(f"\n")

    file.write(f"shell cd dump_file\n"
                "dump 11 all custom 500 dump." + case_name + "_*.cfg mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n\n"
                "reset_timestep 0\n"
                "thermo 10\n"
                "thermo_style custom step pe lx ly lz press\n"
                "min_style cg\n"
                "minimize 1e-25 1e-25 5000 10000\n\n"
                "shell cd ../../\n"
                "shell mkdir result\n"
                "shell cd result\n"
                "dump dump_final all custom 1 "
                "dump_final_" + str(dislocation_line_position[0]*100) + str(dislocation_line_position[1]*100) +
                ".cfg mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
                "run 0\n"
                "undump dump_final\n"
                "write_restart 0K_Emin.restart\n"
                "print"' "'"All done"'"' )

    file.close()

def write_in_file_heat_emin(config, elem_property, directory, shell_radius, sample_center, dislocation_line_position):
  
    pot_path = "/gauss12/home/cityu/anwenliu/git/anisotropic-elasticity/src/anwen_disl_aniso/potential/"
    if elem_property["cryst_structure"] == 'fcc' :
        atom_vol = str(elem_property["latt_const"]**3/4/100000)
    if elem_property["cryst_structure"] == 'hcp' :
        atom_vol = str(np.sqrt(3)*elem_property["latt_const"]["a"]**2*elem_property["latt_const"]["c"]/4/100000)
    elem_symbol = elem_property["elem_symbol"]
    case_name = config["case_name"]
    pot_file_name = elem_property["pot_file_name"]
    mass = str(elem_property["mass"])
 
    filename = 'in.'+ case_name
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"clear\n"
                "units metal\n"
                "dimension 3\n"
                "boundary f f p\n"
                "atom_style atomic\n"
                "atom_modify map array\n\n"
                "read_data " + case_name + "_aniso.dat\n\n"
                "mass 1 " + mass + "\n\n") 
    
    file.write(f"region inner cylinder z {sample_center[0]} {sample_center[1]} {shell_radius} EDGE EDGE\n"
                "group free_atom region inner\n"
                "group freeze_atom subtract all free_atom\n"
                "fix freeze_boundary freeze_atom setforce 0.0 0.0 0.0\n\n"

                "pair_style deepmd frozen_model.pb\n"
                "pair_coeff")
      
    file.write( "\n"
                "neighbor 2.0 bin\n\n"
                "compute freeatomtemp all temp\n"
                "compute peratom all stress/atom NULL\n\n")
  
    for i in range(1,7):
      file.write(f"variable s{lmp_stress_atom[i]} atom c_peratom[{i}]/" + atom_vol + '\n')
    file.write(f"\n")
  
    file.write(f"shell cd dump_file\n"
                "dump 11 all custom 1 dump." + case_name + "_0 mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
                "dump_modify 11 element Ti sort id\n"
                "run 0\n"
                "undump 11\n\n"
                "reset_timestep 0\n"
                "thermo 10\n"
                "thermo_style custom step pe lx ly lz press c_freeatomtemp\n"
                "velocity free_atom create 100.0 4928459 rot yes dist gaussian\n"
		"fix warm free_atom nvt temp 100.0 100.0 0.1\n"
		"run 2000\n"
		"unfix warm\n\n"
		"fix warm free_atom nvt temp 100.0 0.1 0.1\n"
		"run 50000\n"
		"unfix warm\n\n"
                "min_style cg\n"
                "minimize 1e-25 1e-25 5000 10000\n\n"
                "dump 11 all custom 1 dump." + case_name + "_* mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
                "dump_modify 11 element Ti sort id\n"
                "run 0\n"
                "undump 11\n\n"
                "shell cd ../../\n"
                "shell mkdir final_result_compare\n"
                "shell cd final_result_compare\n"
                "dump dump_final all custom 1 "
                "dump.final_" + str(int(dislocation_line_position[0]*100))  + str(int(dislocation_line_position[1]*100)) +
                " mass type x y z v_s11 v_s22 v_s33 v_s12 v_s13 v_s23\n"
                "dump_modify dump_final element Ti sort id\n"
                "run 0\n"
                "undump dump_final\n"
                "print"' "'"All done"'"' )
  
    file.close()



  
def write_job(config, elem_property, directory, core_number):
    '''
    write lammps in file
    '''
    case_name = config["case_name"]
    in_filename = 'in.'+ case_name
    core_num = str(core_number) 
   
    filename = 'job.sh'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')
  
    file.write(f"#!/usr/bin/env bash\n\n"
                "#SBATCH --job-name="'"' + case_name + "_E_min"'"\n'
                "#SBATCH --partition=xlong\n"
                "#SBATCH --ntasks-per-core=1\n"
                "#SBATCH --ntasks=" + core_num + "\n"
                "#SBATCH --nodes=1\n"
                "#SBATCH --exclude=gauss5\n"
                "#SBATCH --cpus-per-task=1\n"
                "#SBATCH --mem=2000M\n"
                "#SBATCH --error='"'j%j.stderr'"'\n"
                "#SBATCH --output='"'j%j.stdout'"'\n\n"
                "# import some basic util func\n"
                ". $HOME/template/shell/bash/bash_script_home\n"
                ". $HOME/template/shell/bash/sys_util\n\n"
                "export OMP_NUM_THREADS=$'{'SLURM_CPUS_PER_TASK'}'\n"
                )
    if elem_property["cryst_structure"] == 'hcp' :
        file.write(f"module load lammps/20200303/aocc_openmpi_znver2_deepmd_xmeam_ipt\n")
    else:
        file.write(f"module load lammps/gcc_openmpi_zen\n")
    file.write(f"ulimit -s unlimited\n\n"
                "appexe='"'nice -n20 lmp'"'\n\n"
                "input_file=" + in_filename +"\n"
                "output_file=$input_file.stdout\n"
                "log_file=$input_file.log\n\n"
                "mpirun -np $SLURM_NTASKS $appexe -in $input_file -log $log_file > $output_file\n")
    file.close()
    
def write_dp_job(config, elem_property, directory, core_number):
    
    case_name = config["case_name"]
    in_filename = 'in.'+ case_name
    core_num = str(core_number)

    filename = 'dpjob.sh'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"#!/usr/bin/env bash\n\n"
                "#SBATCH --job-name="'"' + case_name + "_E_min"'"\n'
                "#SBATCH --partition=gpu_7d1g\n"
                "#SBATCH --nodes=1\n"
                "#SBATCH --ntasks-per-node=1\n"
                "#SBATCH --time=24:00:00\n"
                "#SBATCH --qos=normal\n"
                "#SBATCH --gres=gpu:1\n"
                "#SBATCH --error='"'j%j.stderr'"'\n"
                "#SBATCH --output='"'j%j.stdout'"'\n\n"
                "date\n"
                "#conda activate /home/tongqwen/deepmd-kit-gpu\n"
                "#export OMP_NUM_THREADS=1\n"
                "#export TF_NUM_INTEROP_THREADS=1; export TF_NUM_INTRAOP_THREADS=1;\n"
                "module load openmpi/4.0.5/gcc/8.3.0\n"
                "module load cuda/11.0.2\n"
                "module load cuda/blas/11.0.2\n"
                "module load cuda/fft/11.0.2\n"
                "export LD_LIBRARY_PATH=/cm/shared/apps/cuda11.0/cudnn/lib64:$LD_LIBRARY_PATH\n"
                "#module purge\n"
                "#conda activate\n\n"
                "input_file=in.lammps\n"
                "output_file=$input_file.stdout\n"
                "log_file=$input_file.log\n\n"
                "mpirun -n 1 /home/tongqwen/software/lammps-stable_29Oct2020-gpu/src/lmp_mpi -in $input_file -log $log_file > $output_file\n"
                "date\n")
    file.close()

def write_dp_cpu_job(config, elem_property, directory):
    case_name = config["case_name"]
    in_filename = 'in.'+ case_name

    filename = 'dpjob.sh'
    file_path_name = os.path.join(directory, filename)
    file = open(file_path_name, 'w')

    file.write(f"#!/usr/bin/env bash\n\n"
                "#SBATCH --job-name="'"' + case_name + "_E_min"'"\n'
                "#SBATCH --partition=cpu_14d1n\n"
                "#SBATCH --nodes=1\n"
                "#SBATCH --ntasks-per-core=1\n"
                "#SBATCH --ntasks=128\n"
                "#SBATCH --cpus-per-task=1\n"
                "#SBATCH --mem=480G\n"
                "#SBATCH --error='"'j%j.stderr'"'\n"
                "#SBATCH --output='"'j%j.stdout'"'\n\n"
                "ulimit -s unlimited\n\n"
                "module load lammps/20201029/gcc_intelmpi_deepmd_api_znver2\n"
                "appexe='nice -n20 lmp_mpi'\n\n"
                "input_file=" + in_filename +"\n"
                "output_file=$input_file.stdout\n"
                "log_file=$input_file.log\n\n"
                "mpirun -np $SLURM_NTASKS $appexe -in $input_file -log $log_file > $output_file\n")

    file.close()

