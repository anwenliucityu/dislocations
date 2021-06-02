import numpy as np

def quadrupolar_params(disl_center, box_boundary, repeat_para, unit_cell_size):

    Lx = box_boundary[0][1] - box_boundary[0][0]
    Ly = box_boundary[1][1] - box_boundary[1][0]
    Lz = box_boundary[2][1] - box_boundary[2][0]
    start_x = box_boundary[0][0]
    start_y = box_boundary[1][0]
    start_z = box_boundary[2][0]
    disl_center_1 = [disl_center[0]*Lx+start_x-0.00001,
                     disl_center[1]*Ly+start_y-0.00001]
    if disl_center[0] < 0.5:
        new_center_x = disl_center[0]+0.5
    elif dis_center < 1.0:
        new_center_x = disl_center[0]-0.5
    if disl_center[1] < 0.5:
        new_center_y = disl_center[1]+0.5
    elif dis_center < 1.0:
        new_center_y = disl_center[1]-0.5
    S_area = repeat_para.sum(axis=0)[1] *repeat_para.sum(axis=0)[2]* unit_cell_size[1] *unit_cell_size[0]

    disl_center_2 = [new_center_x*Lx+start_x-0.00001, 
                     new_center_y*Ly+start_y-0.00001]
    d = [0, Ly/2, Lz/2]
    l = [1,0,0]
    A = np.cross(l,d)
    return disl_center_1, disl_center_2, d, l, A
