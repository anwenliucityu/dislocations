'''wrap all atoms into one period by pbc'''

import numpy as np

def pbc_wrap_orthogonal(atom_coor, box_boundary):
    '''
    the box is orthogonal
    '''

    # size of box (Angstrom)
    box_size = [box_boundary[0][1] - box_boundary[0][0],
                box_boundary[1][1] - box_boundary[1][0],
                box_boundary[2][1] - box_boundary[2][0]]

    for axis in range(3):
    	for i in range(1,50):
            tempt = atom_coor[:, axis]
            tempt[np.where(tempt > box_boundary[axis][1])] -= box_size[axis]
            tempt[np.where(tempt < box_boundary[axis][0])] += box_size[axis]
            atom_coor[:, axis] = tempt

    return atom_coor
    
def duplicate_z(atom_coor, atom_type, box_boundary, duplicate_time):
    # duplicate box boundary in z direction
    num_atom = atom_coor.shape[0]
    zlo = box_boundary[2][0]
    zhi = box_boundary[2][1]
    repeat_unit = zhi - zlo
    new_zhi = zlo + repeat_unit * duplicate_time
    new_box_boundary = [[box_boundary[0][0], box_boundary[0][1]],
                        [box_boundary[1][0], box_boundary[1][1]],
                        [zlo, new_zhi]]
    
    new_coor = []
    new_type = []
    shift = np.array([0,0,repeat_unit])
    for i in range(duplicate_time):
        new_coor = np.append(new_coor, np.array(atom_coor) + shift*i)
        new_type = np.append(new_type, atom_type)
    return new_coor.reshape((num_atom*duplicate_time,3)), new_type.reshape((num_atom*duplicate_time,1)), new_box_boundary
