from input_dict import hcp_screw_a_basal as config
from perfect_cryst_constructor import perfect_cryst_constructor, perfect_hcp_constructor
import sample_region
from elast_const_transform import elast_const_transform
from integrand import integrand
from integrate import zero_pi_integrate, zero_theta_integrate
from aniso_disl_theory import initialize_disl_config, stress_field, displacement_field
import pbc_wrap
import numpy as np
from write_output import write_datafile, write_cfg_auxiliary, write_in_file, write_job, output_dir, write_in_finite_T_meta_file
import os
import matplotlib.pyplot as plt
from euler_rotation import atom_torsion
from quadrupolar import quadrupolar_params
import deform_gradient
import multiprocessing

def aniso_disl_constructor(main_path, config_style, dislocation_type,                                                                           
                           pot_element, latt_const, pot_name, element_struct, unit_cell_size,
                           simulation_type, start_temp,
                           elastic_const, mass, pot_cutoff, in_pot, pot_type, pot_path,
                           disl_center,
                           partition, mem, module_load, appexe, ncpu,
                           running_steps=None, dump_interval=None,
                           boundary_freeze_width = 3.0, calc_atomic_stress=True, calc_aniso_stress=False,
                           global_emin=False, temp = 300, spring_factor_k=None, output_perfect=True,
                           cooling_rate=500, sbatch_job=False, ovito=True):

    '''
    main function for constructing an anisotropic dislocation configuration
    '''
    # Creat a new directory path for installing output file
    directory = output_dir(main_path, pot_element, pot_name, element_struct, config_style,
            dislocation_type, simulation_type, unit_cell_size, global_emin=global_emin, T=start_temp, temp=temp)
    
    # construct a perfect crystal
    if element_struct == 'hcp':
        atom_coor, atom_type, box_boundary, new_box_lattice_const, frame_new, repeat_para = \
                perfect_hcp_constructor(config, config_style, latt_const, pot_element, dislocation_type, unit_cell_size, quadru=True)
    if element_struct == 'fcc' or element_struct == 'bcc':
        atom_coor, atom_type, box_boundary = perfect_cryst_constructor(config, config_style,
                element_struct, dislocation_type, latt_const, unit_cell_size,)
        frame_new = config["frame_new"]

    # wrap all atoms into one period by PBC
    atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary) 

    # write perfect configuration for dislocation analysis 
    if output_perfect == True:
        write_datafile(pot_element, dislocation_type, mass, directory, 
            atom_coor, atom_type, box_boundary, suffix='_perfect_ref.dat')
    '''

    # normalize the initial and new coordinate
    frame_initial = config["frame_initial"]
    frame_initial /= np.linalg.norm(frame_initial, axis=0)
    frame_new /= np.linalg.norm(frame_new, axis=0)

    # new elastic constant under new rotated coordinate system
    elastic_const_new = elast_const_transform(elastic_const, frame_initial, frame_new)

    # dislocation parameter
    frame_crystal = config["frame_initial"]
    if element_struct == 'hcp':
        __, b_vector = initialize_disl_config(config, latt_const, element_struct, 
                frame_crystal, new_box_lattice_const, disl_center, repeat_para=repeat_para) 
    else:
        __, b_vector = initialize_disl_config(config, latt_const, element_struct,
                frame_new, latt_const, disl_center) 
    
    # initialize quadrupolar params
    disl_center_1, disl_center_2, d, l, A = quadrupolar_params(disl_center, box_boundary)

    # calculate s_theta, q_theta, S, Q, B, S_theta, Q_theta
    integrate_stepnumber = 100
    theta_list, q, s, b = integrand(elastic_const_new, integrate_stepnumber)
    S, S_list = zero_pi_integrate(s)
    Q, Q_list = zero_pi_integrate(q)
    B, B_list = zero_pi_integrate(b)
    S_theta, s_theta = zero_theta_integrate(s, atom_coor, disl_center_1, S_list)
    Q_theta, q_theta = zero_theta_integrate(q, atom_coor, disl_center_1, Q_list)

    # displacement
    u_displacement_field_1 = displacement_field(atom_coor, disl_center_1, S, B, S_theta, Q_theta, b_vector)
    S_theta, s_theta = zero_theta_integrate(s, atom_coor, disl_center_2, S_list)
    Q_theta, q_theta = zero_theta_integrate(q, atom_coor, disl_center_2, Q_list)
    b_vector[0] = -b_vector[0]
    u_displacement_field_2 = displacement_field(atom_coor, disl_center_2, S, B, S_theta, Q_theta, b_vector) 
    atom_coor = atom_coor + u_displacement_field_1 +u_displacement_field_2
    atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)
  
    # add homogeneous strain to stablized dislocations
    b_vector[0] = -b_vector[0]
    strain = deform_gradient.homo_strain(S_area,b_vector,A)
    deform_grad = deform_gradient.deformation_gradient(strain, F21=0, F31=0, F32=0)
    atom_coor = deform_gradient.deform_atom_coor(deform_grad, atom_coor)
    box_boundary, tilt_para = deform_gradient.deform_box_boundary(deform_grad, box_boundary)
  
    # write datafile
    if True:
        if unit_cell_size[2] == 1:
            write_datafile(pot_element, dislocation_type, mass, directory, 
                atom_coor, atom_type, box_boundary)
        else:
            dup_atom_coor, dup_atom_type, dup_box_boundary = pbc_wrap.duplicate_z(atom_coor,
                    atom_type, box_boundary, unit_cell_size[2]-1)
            write_datafile(pot_element, dislocation_type, mass, directory,
            dup_atom_coor, dup_atom_type, dup_box_boundary)

    # write in file for lammps and return the directory path
    shell_radius = sample_radius - boundary_freeze_width * pot_cutoff # the radius of region of free-moved atoms

    if simulation_type == 'energy_minimization':
        write_in_file(pot_path, config_style, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, shell_radius, sample_center, disl_center, calc_atomic_stress=calc_atomic_stress, global_emin=global_emin,temp=temp,
                  cooling_rate=cooling_rate)
    if simulation_type == 'metastable' or simulation_type == 'finite_T':
        write_in_finite_T_meta_file(pot_path, directory, element_struct, latt_const, pot_element, dislocation_type, in_pot, pot_type,
                  mass, shell_radius, sample_center, disl_center, simulation_type, spring_factor_k=spring_factor_k,
                  running_steps=running_steps,
                  dump_interval=dump_interval, calc_atomic_stress=False, temp = temp)

    # write job.sh file for submitting job
    write_job(pot_element, dislocation_type, partition, mem, module_load, appexe, directory, ncpu=ncpu)

    # ovito
    if ovito == True:
        init_path = os.getcwd()
        os.chdir(directory)
        os.system(f"ovito {pot_element}_{dislocation_type}_disl.dat")
        os.chdir(init_path)

    # cd job.sh directory to submit job
    if sbatch_job == True:
        init_path = os.getcwd()
        os.chdir(directory)
        os.system('sbatch job.sh')
        os.chdir(init_path)
    '''
    return directory

if __name__ == '__main__':
    aniso_disl_constructor()
