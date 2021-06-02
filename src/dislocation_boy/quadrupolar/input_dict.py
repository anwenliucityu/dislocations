'''input file'''

import math

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
    "frame_new"         : [[   1,  1.0, -0.5],
                           [  -1,  1.0, -0.5],
                           [   0,  1.0,  1.0]],
}

Al_screw_config = {
    "case_name"         : "Al_screw",
    "elem"              : "EAM_Mishin_Al",
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
    "type_to_elem"      : {1 : "EAM_Mishin_Al", 2 : "EAM_Mishin_Al"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 100,#52,
    "box_num_cell_y"    : 75,#36,
    "box_num_cell_z"    : 50,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [-0.5, 0.5, 0.0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [-1./math.sqrt(2.), 1./math.sqrt(2.), 0.],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.5,
    "disl_center_y"     : 0.5,
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

Ti_screw_basal_a_config = {
    "case_name"         : "Ti_screw_basal_a",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 6,
    "box_num_cell_y"    : 8,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0., 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   0,  -1 ],
                           [   0,   1,  0 ]],
}

Ti_edge_basal_a_config = {
    "case_name"         : "Ti_edge_basal_a",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 180,
    "box_num_cell_y"    : 120,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   0,  -1 ],
                           [   0,   1,  0 ]],
}

Ti_screw_prismI_a_config = {
    "case_name"         : "Ti_screw_prismI_a",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 130,
    "box_num_cell_y"    : 120,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   0,  0,  -1 ],
                           [   0,  1,  0 ],
                           [   1,  0,  0 ]],
}

Ti_edge_prismI_a_config = {
    "case_name"         : "Ti_edge_prismI_a",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 175,
    "box_num_cell_y"    : 115,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1, 0, 0],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
    # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   0,  -1 ],
                           [   0,   1,  0 ]],
}

Ti_edge_pyrI_a_config = {
    "case_name"         : "Ti_edge_pyrI_a",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 200,
    "box_num_cell_y"    : 140,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    "frame_new"         :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],
}

Ti_screw_pyrII_ac_config = {
    "case_name"         : "Ti_screw_pyrII_ac",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 100,
    "box_num_cell_y"    : 200,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
     # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],

}

Ti_edge_pyrII_ac_config = {
    "case_name"         : "Ti_edge_pyrII_ac",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 110,
    "box_num_cell_y"    : 240,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],

}

Ti_mixed_prismI_ac_config = { #NO PBC
    "case_name"         : "Ti_mixed_prismI_ac",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 180,
    "box_num_cell_y"    : 120,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],

} # no pbc !!>>>

Ti_mixed_pyrI_ac_config = { #NO PBC
    "case_name"         : "Ti_mixed_pyrI_ac",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 50,
    "box_num_cell_y"    : 120,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 0, -1],
                            [0, 1, 0]],
}


Ti_screw_prismI_c_config = {
    "case_name"         : "Ti_screw_prismI_c",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 175,
    "box_num_cell_y"    : 105,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [0, 0,  1],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],

}

Ti_edge_prismI_c_config = {
    "case_name"         : "Ti_edge_prismI_c",
    "elem"              : "XMEAM_Rui_Ti",
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
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 130,
    "box_num_cell_y"    : 125,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[0, 0, -1],
                            [0, 1, 0],
                            [1, 0, 0]],

}

Ti_edge_prismII_c_config = {
    "case_name"         : "Ti_edge_prismII_c",
    "elem"              : "XMEAM_Rui_Ti",
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
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 120,
    "box_num_cell_y"    : 190,
    "box_num_cell_z"    : 1,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1, 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [0, 0, 1],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./12,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
     "frame_new"        :  [[0, 0, -1],
                            [0, 1, 0],
                            [1, 0, 0]],

}

Ti_screw_basal_a_quadrupole_config = {
    "case_name"         : "quadru_Ti_screw_basal_a",
    "elem"              : "XMEAM_Rui_Ti",
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
    "type_basis_atoms"  : [1, 1, 1, 1,],
    # bind type to element
    "type_to_elem"      : {1 : "XMEAM_Rui_Ti", 2 : "XMEAM_Rui_Ti"},
    # number of replicate of unit cell
    "box_num_cell_x"    : 1,
    "box_num_cell_y"    : 6,
    "box_num_cell_z"    : 8,
    # Burgers vector (scaled by lattice constant)
    "burgers"           : [1., 0,  0],
    # dislocation line direction (unit vector)
    "disl_line_direction" : [1., 0, 0],
    # dislocation position (relative coordinates in a unit cell)
    "disl_center_x"     : 0.00001,
    "disl_center_y"     : -1./6,
        # initial frame: [e1, e2, e3], ei are column vectors
    "frame_initial"     :  [[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]],
    # new frame: [e1', e2', e3'], ei' are column vectors
    "frame_new"         : [[   1,   0,  0 ],
                           [   0,   1,  0 ],
                           [   0,   0,  1 ]],
}
