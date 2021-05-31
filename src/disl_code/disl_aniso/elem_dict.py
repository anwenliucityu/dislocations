'''list of properties of elements'''

elem_dict = {
    "EAM_Mishin_Al" : {
        "elem_symbol"       : 'Al',
        "pot_file_name"     : ['Al99.eam.alloy'],
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 4.05000466178543,      # lattice constant (Angstrom)
        "mass"              : 26.98,
        "pot_cutoff"        : 6.287,     # cut-off radius of potential (Angstrom)
        "C11"               : 1.14e11,   # (Pa)
        "C12"               : 0.619e11,
        "C44"               : 0.316e11,
        "shear_modulus"     : 0.316e11,  # C44
        "poisson_ratio"     : 0.331016,  # C12/(2 * (C44 + C12))                                                                           
        "elastic_const_initial"  :  [1.14e11, 0.619e11, 0.619e11,   0,   0,   0,
                                               1.14e11, 0.619e11,   0,   0,   0,
                                                         1.14e11,   0,   0,   0,
                                                             0.316e11,   0,   0,
                                                                  0.316e11,   0,
                                                                       0.316e11],
    },

    "EAM1_Mishin_Cu" : {
        "elem_symbol"       : 'Cu',
        "pot_file_name"     : ['Cu01.eam.alloy'],
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 3.615,     # lattice constant (Angstrom)
        "mass"              : 63.546,
        "pot_cutoff"        : 5.50679,   # cut-off radius of potential (Angstrom)
        "C11"               : 1.699e11,  # (Pa)
        "C12"               : 1.226e11,
        "C44"               : 0.762e11,
        "shear_modulus"     : 0.762e11,  # C44
        "poisson_ratio"     : 0.3083501, # C12/(2 * (C44 + C12))
        "elastic_const_initial"  :  [1.699e11, 1.226e11, 1.226e11,   0,   0,   0,
                                               1.699e11, 1.226e11,   0,   0,   0,
                                                         1.699e11,   0,   0,   0,
                                                              0.762e11,   0,   0,
                                                                   0.762e11,   0,
                                                                        0.762e11],

    },
        "EAM_Mishin_Ni" : {
        "elem_symbol"       : 'Ni',
        "pot_file_name"     : ['Ni99.eam.alloy'],
        "pair_style"        : 'eam/alloy',
        "cryst_structure"   : 'fcc',     # crystal structure
        "latt_const"        : 3.52,     # lattice constant (Angstrom)
        "mass"              : 58.71,
        "pot_cutoff"        : 5.804,   # cut-off radius of potential (Angstrom)
        "C11"               : 2.47e11,  # (Pa)
        "C12"               : 1.48e11,
        "C44"               : 1.25e11,
        "shear_modulus"     : 1.25e11,  # C44
        "poisson_ratio"     : 0.27106227106227104, # C12/(2 * (C44 + C12))
        "elastic_const_initial"  :  [2.47e11,  1.48e11,  1.48e11,   0,   0,   0,
                                               2.47e11,  1.48e11,   0,   0,   0,
                                                         2.47e11,   0,   0,   0,
                                                              1.25e11,   0,   0,
                                                                   1.25e11,   0,
                                                                        1.25e11],

    },
        "Rui_Ti" : {
        "elem_symbol"       : 'Ti',
        "pot_file_name"     : ['xmeam_Ti/library.xmeam', 'xmeam_Ti/Ti.xmeam'],
        "pair_style"        : 'xmeam',
        "cryst_structure"   : 'hcp',     # crystal structure
        "latt_const"        : { "a" : 2.95436648798772, "c" : 2.95436648798772*1.58282966743688},    # lattice constant a & c (Angstrom)
        "mass"              : 47.867 ,
        "pot_cutoff"        : 7.4602,   # cut-off radius of potential (Angstrom)
        "elastic_const_initial"  :  [1.80576472663722e11,  0.916713001972235e11,  0.671350317864564e11,   2.63112173751293e-09,   1.50792656841712e-10,   -0.00000400735785371715,
                                               1.8057647558547e11,  0.67135030672858e11,   -2.64119117391231e-09,   -1.5000676442797e-10,   0.0000040115416852069,
                                                         1.9267750723697e11,   -6.99294451928992e-12,   -1.42852448965038e-11,   -1.4611200551433e-09,
                                                              0.322381231567968e11,   -5.26135002680259e-10,   -2.9673473153242e-10,
                                                                   0.322381235556176e11,   1.47357562711966e-09,
                                                                        0.444512187797651e11],
     },
        "XMEAM_Rui_Ti" : {#DP potential
        "elem_symbol"       : 'Ti',
        "pot_file_name"     : ['xmeam_Ti/library.xmeam', 'xmeam_Ti/Ti.xmeam'],
        "pair_style"        : 'xmeam',
        "cryst_structure"   : 'hcp',     # crystal structure
        "latt_const"        : { "a" : 2.93849, "c" :4.65842 },    # lattice constant a & c (Angstrom) at 300K
        "mass"              : 47.867 ,
        "pot_cutoff"        : 7.4602,   # cut-off radius of potential (Angstrom)
        "elastic_const_initial"  :  [1.74483121359734e11,  0.799737260851811e11,  0.820045532917221e11,   3.7344994750608e-15,   -1.02942178102389e-15,   7.03084821743516e-15,
                                                           1.74333257466625e11,   0.819688541344252e11,   -1.44897141527536e-14,   7.58600937367035e-15,   -3.40938063011528e-14,
                                                                                  2.00359366354643e11,   -7.32666799108776e-15,   -1.54586491458141e-15,   5.75029746180558e-14,
                                                                                                          0.438751819396401e11,   6.149369629665e-18,   -6.12852430888515e-17,
                                                                                                                                    0.43875189892828e11,   2.56166062515217e-14,
                                                                                                                                                           0.472792176422276e11],

    },
        "MEAM_Wu_Mg" :{
        "elem_symbol"       : 'Mg',
        "pot_file_name"     : ['meam_Mg/library.meam', 'meam_Mg/Mg.meam'], ## 0K
        "pair_style"        : 'meam/c',
        "cryst_structure"   : 'hcp', 
        "latt_const"        : { "a" : 3.18553836632149, "c" :5.16735066638725 },
        "mass"              : 24.305,
        "pot_cutoff"        : 8.2901, 
        "elastic_const_initial"  : [64.2141858981347e9, 26.6954735937704e9, 21.9647112734493e9, 0,0,0,
                                                        64.2141858932623e9, 21.9647112713084e9, 0,0,0,
                                                                            70.0829987796629e9, 0,0,0,
                                                                               17.7656370295554e9,0,0,
                                                                                 17.7656370293065e9,0,
                                                                                   18.7593559578235e9]
        
        },
        
        "MEAM_Wu_Mg_300K" :{ #300K
        "elem_symbol"       : 'Mg',
        "pot_file_name"     : ['meam_Mg/library.meam', 'meam_Mg/Mg.meam'], 
        "pair_style"        : 'meam/c',
        "cryst_structure"   : 'hcp', 
        "latt_const"        : { "a" : 3.2056898354662, "c" :5.20332683029676 },
        "mass"              : 24.305,
        "pot_cutoff"        : 8.2901, 
        "elastic_const_initial"  : [57.4970073278751e9, 25.3908076752674e9, 19.9233473899611e9, 0,0,0,
                                                        56.9582458612214e9, 20.1979299347709e9, 0,0,0,
                                                                            63.398365797924e9, 0,0,0,
                                                                               14.7539188100717e9,0,0,
                                                                                 15.1426107230409e9,0,
                                                                                   16.0219401091426e9]

        
        },
        "XMEAM_Rui_Mg":{
        "elem_symbol"       : 'Mg',
        "pot_file_name"     : ['xmeam_Mg/library.xmeam', 'xmeam_Mg/Mg.xmeam'], ## 0K
        "pair_style"        : 'xmeam',
        "cryst_structure"   : 'hcp', 
        "latt_const"        : { "a" : 3.2056898354662, "c" :5.20332683029676 },
        "mass"              : 24.305,
        "pot_cutoff"        : 6.1, 
        "elastic_const_initial"  : [64.3022564913879e9, 25.3591979096868e9, 21.7643309570695e9, 0,0,0,
                                                        64.3054156624858e9, 21.7642366072362e9, 0,0,0,
                                                                            69.844854376388e9, 0,0,0,
                                                                               17.977698852538e9,0,0,
                                                                                 17.9779306271668e9,0,
                                                                                   19.4772889693936e9]
        
        },
        "MEAM_Hennig_bcc_Ti":{ #0K
        "elem_symbol"       : 'Ti',
        "pot_file_name"     : ['Ti.meam.spline'],
        "pair_style"        :'meam/spline',
        "cryst_structure" : 'bcc',
        "latt_const"        : 3.27156976391928,
        "mass"              : 47.867,
        "pot_cutoff"        : 8,
        
        "elastic_const_initial"  : [94.8461044876344e9, 110.973214688863e9, 110.973214688863e9, 0,0,0,
                                                        94.8461044876344e9, 110.973214688863e9, 0,0,0,
                                                                            94.8461044876344e9, 0,0,0,
                                                                               52.8631997890871e9,0,0,
                                                                                 52.8631997890871e9,0,
                                                                                   52.8631997890871e9]


        }
}
