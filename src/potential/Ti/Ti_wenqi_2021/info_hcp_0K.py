#!/usr/bin/python3.8 

#########################################
#      information of DP potential      #
#########################################

# information of material
'''
    units:
        temperature      : Kelvin
        lengh            : Angstrom
        mass             : g/cm^3
        elastic constant : Pa
'''
latt_const     =    {'a': 2.93749682368896, 'c': 4.6463868242434}
mass           =    47.867
c11            =    1.74483121359734e11
c22            =    1.74333257466625e11
c33            =    2.00359366354643e11
c12            =    0.799737260851811e11
c13            =    0.820045532917221e11
c23            =    0.819688541344252e11
c44            =    0.438751819396401e11
c55            =    0.43875189892828e11
c66            =    0.472792176422276e11
c14            =    3.7344994750608e-6
c15            =    -1.02942178102389e-6
c16            =    7.03084821743516e-6
c24            =    -1.44897141527536e-5
c25            =    7.58600937367035e-6
c26            =    -3.40938063011528e-5
c34            =    -7.32666799108776e-6
c35            =    -1.54586491458141e-6
c36            =    5.75029746180558e-5
c45            =    6.149369629665e-9
c46            =    -6.12852430888515e-8
c56            =    2.56166062515217e-5

# initialize elastic constant
elastic_const  = [c11, c12, c13, c14, c15, c16,
                       c22, c23, c24, c25, c26,
                            c33, c34, c35, c36,
                                 c44, c45, c46,
                                      c55, c56,
                                           c66]




