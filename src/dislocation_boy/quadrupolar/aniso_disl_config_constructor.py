from input_dict import Ti_screw_basal_a_config as config
from perfect_cryst_constructor import perfect_cryst_constructor, perfect_hcp_constructor
import sample_region
from elast_const_transform import elast_const_transform
from integrand import integrand
from integrate import zero_pi_integrate, zero_theta_integrate
from aniso_disl_theory import initialize_disl_config, stress_field, displacement_field
import pbc_wrap
from elem_dict import elem_dict
import numpy as np
from write_output import write_datafile, write_cfg_auxiliary, write_dp_in_file, write_dp_job, write_info, output_directory, write_in_file_heat_emin
import write_output as wo
import os
import time
from check_lmp import check_lmp
from euler_rotation import atom_torsion
import multiprocessing
from euler_rotation import Rx, Ry, Rz
import deform_gradient

def aniso_disl_constructor():
  '''
  main function for constructing an anisotropic dislocation configuration
  '''
  #-------------------
  start = time.time()
  #-------------------

  # material property
  elem_property = elem_dict[config["elem"]]
 
  # Creat a new directory path for writing output file
  project_name = str("Ti_DFT_model_oppo_burger_accurate_dis_line")
  path = '/gauss12/home/cityu/anwenliu/scratch/dislocation_simulation/' + project_name
  directory = output_directory(config, elem_property, path)
   
  # construct a perfect crystal
  atom_coor, atom_type, box_boundary, lattice_const, frame_new, repeat_para, tilt, S_area = perfect_hcp_constructor(config, elem_property) 
  # wrap all atoms into one period by PBC
  atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)
    
  # normalize the initial and new coordinate
  frame_initial = config["frame_initial"]
  frame_initial /= np.linalg.norm(frame_initial, axis=0)
  frame_new /= np.linalg.norm(frame_new, axis=0)

  # new elastic constant under new rotated coordinate system
  elastic_const_initial = np.array(elem_property["elastic_const_initial"])
  elastic_const_new = elast_const_transform(elastic_const_initial, frame_initial, frame_new)
  
  # dislocation parameter
  frame_crystal = config["frame_initial"]
  disl_center, b_vector = initialize_disl_config(config, elem_property, frame_crystal, lattice_const) # frame_crystal changed to frame new for fcc
  Lx = box_boundary[0][1] - box_boundary[0][0]
  Ly = box_boundary[1][1] - box_boundary[1][0]
  start_x = box_boundary[0][0]
  start_y = box_boundary[1][0]
  a = 2.9374968236889
  c = 4.6463868242434
  #### randam
  import random
  
  randx = random.random()-0.5
  a1 = np.sqrt(3)*a/2*randx
  randy = random.random()-0.5
  a2 = randy*c
  print(a1,a2)
  
  #disl_center_1 = [0.25*Lx+start_x+np.sqrt(3)*a/6, 0.25*Ly+start_y+c/4]
  #disl_center_2 = [0.75*Lx+start_x+np.sqrt(3)*a/6, 0.75*Ly+start_y+c/4]
  #disl_center_1 = [1.5*np.sqrt(3)*a+1/6*np.sqrt(3)*a+start_x , 2.25*c+start_y] ### c and d site
  #disl_center_2 = [1.5*np.sqrt(3)*a+1/6*np.sqrt(3)*a + 0.5*Lx+start_x , 1.75*c + 0.5*Ly+start_y]
  #disl_center_1 = [1.75*np.sqrt(3)*a+start_x+0.001, 2*c+start_y+0.001]  ## E config
  #disl_center_2 = [1.75*np.sqrt(3)*a+start_x+ 0.5*Lx+0.001, 2*c+start_y+ 0.5*Ly+0.001]
  disl_center_1 = [1.5*np.sqrt(3)*a + a1 + start_x, 2*c +a2+ start_y] ##b config
  disl_center_2 = [1.5*np.sqrt(3)*a + a1 + 0.5*Lx + start_x, 2*c + a2+ start_y + 0.5*Ly]
  print('center = ', disl_center_1)
  #disl_center_1 = [0.25*Lx+start_x+0.01, 0.25*Ly+start_y+0.01]
  #disl_center_2 = [0.75*Lx+start_x+0.01, 0.75*Ly+start_y+0.01]

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
  b_vector[2] = -b_vector[2]
  u_displacement_field_2 = displacement_field(atom_coor, disl_center_2, S, B, S_theta, Q_theta, b_vector)
  atom_coor = atom_coor + u_displacement_field_1 +u_displacement_field_2
  
  # wrap all atoms into one period by PBC
  atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)
  
  # Rotate the coodinate
  
  r_angle = np.pi/2
  rotation_matrix = np.dot(Rx(r_angle),Ry(r_angle))
  
  atom_coor = deform_gradient.rotate_atom_coor(atom_coor, rotation_matrix)
  box_boundary = deform_gradient.rotate_box_boundary(box_boundary, rotation_matrix)
  
  # add homogeneous strain to stablized dislocations
  Lx = box_boundary[0][1] - box_boundary[0][0]
  Ly = box_boundary[1][1] - box_boundary[1][0]
  Lz = box_boundary[2][1] - box_boundary[2][0]
  d = [0, Ly/2, Lz/2]
  l = [1,0,0]
  A = np.cross(l,d)
  b_vector = [-b_vector[2], b_vector[0], b_vector[1]]
  #b_vector = [-b_vector[2], b_vector[0], b_vector[1]]
  strain = deform_gradient.homo_strain(S_area,b_vector,A)
  print(strain)
  
  deform_grad = deform_gradient.deformation_gradient(strain, F21=0, F31=0, F32=0)
  print(deform_grad)
  
  atom_coor = deform_gradient.deform_atom_coor(deform_grad, atom_coor)
  box_boundary, tilt_para = deform_gradient.deform_box_boundary(deform_grad, box_boundary)
  
  # write datafile
  write_datafile(config, directory, atom_coor, atom_type, box_boundary,tilt_para)
  
  # construct a perfect crystal
  atom_coor, atom_type, box_boundary, lattice_const, frame_new, repeat_para, tilt, S_area = perfect_hcp_constructor(config, elem_property) 
  # wrap all atoms into one period by PBC
  atom_coor = pbc_wrap.pbc_wrap_orthogonal(atom_coor, box_boundary)
  
  # Rotate the coodinate                                                                                        
  r_angle = np.pi/2
  rotation_matrix = np.dot(Rx(r_angle),Ry(r_angle))
  atom_coor_ref = deform_gradient.rotate_atom_coor(atom_coor, rotation_matrix)
  box_boundary_ref = deform_gradient.rotate_box_boundary(box_boundary, rotation_matrix)
  
  # write perfect configuration for dislocation analysis 
  wo.write_perfect_datafile(config, directory, atom_coor_ref, atom_type, box_boundary_ref)

  # write in file for lammps and return the directory path
  sample_radius = 100
  shell_radius = sample_radius - 3 * elem_property["pot_cutoff"] # the radius of region of free-moved atoms
  
  #write_in_file_heat_emin(config, elem_property, directory, shell_radius, sample_center, dislocation_line_position)
  #write_dp_in_file(config, elem_property, directory, shell_radius, sample_center, dislocation_line_position)
  sample_center = [0,0]
  dislocation_line_position = [0,0]
  #write_in_file_heat_emin(config, elem_property, directory, shell_radius, sample_center, dislocation_line_position)
  
  # write job.sh file for submitting job
  #wo.write_dp_cpu_job(config, elem_property, directory) # core number
  
  return directory

if __name__ == '__main__':
  aniso_disl_constructor()
