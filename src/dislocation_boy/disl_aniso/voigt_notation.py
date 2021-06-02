'''Voigt notation'''

voigt_one_to_two = {
    0      :   (0, 0),
    1      :   (1, 1),
    2      :   (2, 2),
    3      :   (1, 2),
    4      :   (2, 0),
    5      :   (0, 1),
}

voigt_two_to_one = {
    (0, 0) :    0,
    (1, 1) :    1,
    (2, 2) :    2,
    (1, 2) :    3,      (2, 1) :    3,
    (2, 0) :    4,      (0, 2) :    4,
    (0, 1) :    5,      (1, 0) :    5,
}

four_tensor_two_to_one = {
    (0, 0) :    0,
    (0, 1) :    1,      (1, 0) :    1,
    (0, 2) :    2,      (2, 0) :    2,
    (0, 3) :    3,      (3, 0) :    3,
    (0, 4) :    4,      (4, 0) :    4,
    (0, 5) :    5,      (5, 0) :    5,
    (1, 1) :    6,
    (1, 2) :    7,      (2, 1) :    7,
    (1, 3) :    8,      (3, 1) :    8,
    (1, 4) :    9,      (4, 1) :    9,
    (1, 5) :    10,     (5, 1) :    10,
    (2, 2) :    11,
    (2, 3) :    12,     (3, 2) :    12,
    (2, 4) :    13,     (4, 2) :    13,
    (2, 5) :    14,     (5, 2) :    14,
    (3, 3) :    15,
    (3, 4) :    16,     (4, 3) :    16,
    (3, 5) :    17,     (5, 3) :    17,
    (4, 4) :    18,
    (4, 5) :    19,     (5, 4) :    19,
    (5, 5) :    20,
}

four_tensor_one_to_two = {
    0   :   (0, 0),
    1   :   (0, 1),
    2   :   (0, 2),
    3   :   (0, 3),
    4   :   (0, 4),
    5   :   (0, 5),
    6   :   (1, 1),
    7   :   (1, 2),
    8   :   (1, 3),
    9   :   (1, 4),
    10  :   (1, 5),
    11  :   (2, 2),
    12  :   (2, 3),
    13  :   (2, 4),
    14  :   (2, 5),
    15  :   (3, 3),
    16  :   (3, 4),
    17  :   (3, 5),
    18  :   (4, 4),
    19  :   (4, 5),
    20  :   (5, 5),
}

four_tensor_list_to_matrix = {
    0   :   (0, 0),
    1   :   (0, 1),
    2   :   (0, 2),
    3   :   (1, 0),
    4   :   (1, 1),
    5   :   (1, 2),
    6   :   (2, 0),
    7   :   (2, 1),
    8   :   (2, 2),
}

four_tensor_matrix_to_list = {
    (0, 0)   :   0,
    (0, 1)   :   1,
    (0, 2)   :   2,
    (1, 0)   :   3,
    (1, 1)   :   4,
    (1, 2)   :   5,
    (2, 0)   :   6,
    (2, 1)   :   7,
    (2, 2)   :   8,
}

lmp_stress_atom = {
     1     :     11,
     2     :     22,
     3     :     33,
     4     :     12,
     5     :     13,
     6     :     23,
}