'''input file'''

import math


bcc_Ti_screw_config = {
    "case_name"         : "bcc_Ti_screw",
    "elem"              : "MEAM_Hennig_bcc_Ti",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : math.sqrt(2), #[-1,1,0]
    "cell_y"            : math.sqrt(6)/3,      #[1,1,2]
    "cell_z"            : math.sqrt(3)/2,#[1,1,1]
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [0.5, 0.5, 1./3 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "MEAM_Hennig_bcc_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 40,#44,
    "box_num_cell_y"    : 70,#18,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0.5, 0.5, 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(3.), 1./math.sqrt(3.), 1./math.sqrt(3.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[  -1.0,  -0.5,  1.0],
                           [   1.0,  -0.5,  1.0],
                           [   0,   1.0,  1.0]],
}

Al_edge_config = {
    "case_name"         : "Al_edge",
    "elem"              : "EAM_Mishin_Al",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 140,#44,
    "box_num_cell_y"    : 90,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(6.), -1./math.sqrt(6.), 2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[  -1.0,  1.0, -0.5],
                           [   1.0,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

fcc_screw = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(6),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(2),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [2./3.,  2./3.,  0      ],
                           [1./6.,  2./3.,  1./2.  ],
                           [1./3.,  1./3.,  0      ],
                           [5./6.,  1./3.,  1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(2.), 1./math.sqrt(2.), 0.],
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[ 0.5,  1.0,  1.0],
                           [ 0.5,  1.0, -1.0],
                           [-1.0,  1.0,  0.0]],
}

Al_mix_config = {
    "case_name"         : "Al_mix",
    "elem"              : "EAM_Mishin_Al",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 88,#44,
    "box_num_cell_y"    : 36,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0., 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Cu_edge_config = {
    "case_name"         : "Cu_edge",
    "elem"              : "EAM1_Mishin_Cu",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM1_Mishin_Cu", 2 : "EAM_Mishin_Cu"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 140,#44,
    "box_num_cell_y"    : 58,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0.5, -0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Cu_screw_config = {
    "case_name"         : "Cu_screw",
    "elem"              : "EAM1_Mishin_Cu",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(6),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(2),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [2./3.,  2./3.,  0      ],
                           [1./6.,  2./3.,  1./2.  ],
                           [1./3.,  1./3.,  0      ],
                           [5./6.,  1./3.,  1./2.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM1_Mishin_Cu", 2 : "EAM1_Mishin_Cu"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 103,#52,
    "box_num_cell_y"    : 72,#36,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(2.), 1./math.sqrt(2.), 0.],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[ 0.5,  1.0,  1.0],
                           [ 0.5,  1.0, -1.0],
                           [-1.0,  1.0,  0.0]],
}

Cu_mix_config = {
    "case_name"         : "Cu_mix",
    "elem"              : "EAM1_Mishin_Cu",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM1_Mishin_Cu", 2 : "EAM1_Mishin_Cu"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 176,#44,
    "box_num_cell_y"    : 80,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0., 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Ni_mix_config = {
    "case_name"         : "Ni_mix",
    "elem"              : "EAM_Mishin_Ni",
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 0.5 * math.sqrt(2),
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 0.5 * math.sqrt(6),
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,      0      ],
                           [1./2.,  0,      1./2.  ],
                           [0,      2./3.,  1./3.  ],
                           [1./2.,  2./3.,  5./6.  ],
                           [0,      1./3.,  2./3.  ],
                           [1./2.,  1./3.,  1./6.  ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1],
    # bind type to element
    "type_to_elem"      : {1 : "EAM_Mishin_Ni", 2 : "EAM_Mishin_Ni"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 198,#44,
    "box_num_cell_y"    : 81,#18,
    "box_num_cell_z"    : 2,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0., 0.5],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1./math.sqrt(6.), 1./math.sqrt(6.), -2./math.sqrt(6.)],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

hcp_screw_a_basal = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : math.sqrt(3),
    "cell_y"            : 1,
    "cell_z"            : 1,
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,      0,     0    ],
                           [1./2.,  0,     1./2 ],
                           [1./3,   1./2,  0.   ],
                           [5./6,   1./2,  1./2 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1,],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0., 0, 1],
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   0,  -1 ],
                           [   0,   1,  0 ]],
}

hcp_edge_a_basal = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : math.sqrt(3),
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,         0,     0    ],
                           [1./2.,     0,     1./2 ],
                           [0,      1./2,     2./3 ],
                           [1./2,   1./2,  1./6 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1,],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   0,  -1 ],
                           [   0,   1,  0 ]],
}

hcp_screw_a_prismI = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 1,
    "cell_x_latt_const" : "c",
    "cell_y_latt_const" : "a",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,         0,     0    ],
                           [0,      1./2,     1./2 ],
                           [1./2,   1./3,     0    ],
                           [1./2,   5./6,     1./2 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1,],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
     # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   0,  0,  -1 ],
                           [   0,  1,  0 ],
                           [   1,  0,  0 ]],
}

hcp_edge_a_prismI = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : math.sqrt(3),
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,         0,     0    ],
                           [1./2,      0,     1./2 ],
                           [0,      1./2,     1./3 ],
                           [1./2,   1./2,     5./6 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1,],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1, 0, 0],
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   0,  -1 ],
                           [   0,   1,  0 ]],
}

hcp_edge_a_pyrI = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 2,
    "cell_z"            : math.sqrt(3),
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[1.,        -1.,     1    ],
                           [1./2,     -1.,     3./2 ],
                           [1.,      -3./4,     4./3 ],
                           [1./2,   -3./4,     11./6 ],
                           [1.,      -1./2,     1.    ],
                           [1./2,   -1./2,     1./2 ],
                           [1./2,   -1./4,    5./6 ],
                           [1.,      -1./4,    1./3 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1, 1, 1],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    "frame_new"         :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],
}

hcp_screw_ca_pyrII = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : math.sqrt(3),
    "cell_y"            : 1,
    "cell_z"            : 1,
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,        -1,     1 ],
                           [1./2,     -1,  3./2 ],
                           [1./3,  -1./2,     0 ],
                           [5./6,  -1./2,  1./2 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
     # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],

}

hcp_edge_ca_pyrII = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : math.sqrt(3),
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,          0,     0 ],
                           [-1./2,     0,    1./2 ],
                           [ 0 ,      1./2,  1./3 ],
                           [ 1./2 ,   1./2,  5./6 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],
}

hcp_mixed_ca_prismI = { #NO PBC
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 1,
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "a",
    "cell_z_latt_const" : "c",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,          0,     0 ],
                           [-1./2,    1./2,    0 ],
                           [ 1./2 ,   1./6,  1./2 ],
                           [ 0 ,      2./3,  1./2 ],],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1,],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],

} # no pbc !!>>>

hcp_mixed_ca_pyrI = { #NO PBC
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : math.sqrt(3),
    "cell_y"            : 2,
    "cell_z"            : 1,
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "c",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,          0,        0 ],
                           [ 1./2,      0,     1./2 ],
                           [ 1./3 ,   1./4,       0 ],
                           [ 5./6 ,   1./4,    1./2 ],
                           [ -1./2,   1./2,    1./2 ],
                           [-1./6,    3./4,    1./2 ],
                           [-2./3,    3./4,    0    ],
                           [0,        1./2,    0]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1, 1, 1, 1, 1, ],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],
}

hcp_screw_c_prismI = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 1,
    "cell_x_latt_const" : "a",
    "cell_y_latt_const" : "a",
    "cell_z_latt_const" : "c",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,        0,     0 ],
                           [1./2,   1./2,    0 ],
                           [1./2,   1./6,  1./2 ],
                           [1.,     2./3,  1./2 ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
}

hcp_edge_c_prismI = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : math.sqrt(3),
    "cell_z"            : 1,
    "cell_x_latt_const" : "c",
    "cell_y_latt_const" : "a",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,        0,     0 ],
                           [0,   1./2,    1./2 ],
                           [1./2,   5./6,  1./2 ],
                           [1./2,     1./3,  1. ]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1],
    # bind type to element
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[0, 0, -1],
                            [0, 1, 0],
                            [1, 0, 0]],
}

hcp_edge_c_prismII = {
    # edge length of a unit cell (scaled by lattice constant)
    "cell_x"            : 1,
    "cell_y"            : 1,
    "cell_z"            : math.sqrt(3),
    "cell_x_latt_const" : "c",
    "cell_y_latt_const" : "a",
    "cell_z_latt_const" : "a",
    # basis atoms (relative coordinates)
    "basis_atoms"       : [[0,        0,     0 ],
                           [0,   1./2,    1./2 ],
                           [1./2,  1./2,   5./6 ],
                           [1./2,  1. ,    1./3]],
    # type of each atom
    "type_basis_atoms"  : [1, 1, 1, 1],
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[0, 0, -1],
                            [0, 1, 0],
                            [1, 0, 0]],
}

