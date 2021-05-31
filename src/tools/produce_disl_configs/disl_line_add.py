import sys
sys.path.append('..')
import importlib

# open main function file

config_list = ['Ti_screw_basal_a_config',  'Ti_edge_basal_a_config',
               'Ti_screw_prismI_a_config', 'Ti_edge_prismI_a_config',
               'Ti_edge_pyrI_a_config',    'Ti_screw_pyrII_ac_config',
               'Ti_edge_pyrII_ac_config',  'Ti_screw_prismI_c_config',
               'Ti_edge_prismI_c_config',  'Ti_edge_prismII_c_config',]

dislocation_line_position = [(1./4, 1./4),
                             (3./4, 1./4),
                             (1./4, 3./4),
                             (3./4, 3./4)]

dirs=[]
for i in config_list:
    # open file and change config
    with open("/gauss12/home/cityu/anwenliu/git/anisotropic-elasticity/src/anwen_disl_aniso/aniso_disl_config_constructor.py", 'r') as f1:
        lines = f1.readlines()
    line_zero = lines.pop(0)
    line_zero = line_zero.replace(line_zero.split()[3], i)

    # write new file with new config
    with open("/gauss12/home/cityu/anwenliu/git/anisotropic-elasticity/src/anwen_disl_aniso/aniso_disl_config_constructor.py", 'w') as f2:
        f2.write(line_zero)
        for line in lines:
            f2.write(line)
    # for each config, pick up 4 dependent dislocation line and get the config
    for j in dislocation_line_position:
        import aniso_disl_config_constructor
        importlib.reload(aniso_disl_config_constructor)
        directory = aniso_disl_config_constructor.aniso_disl_constructor(j)
        directory = '\"' + directory + '\"'
        dirs.append(directory)
'''
with open("/gauss12/home/cityu/anwenliu/dislocation/aniso_output/heat_up_emin/directory_list_sbatch.py", 'w') as f3:
    dir_list = ",".join(dirs)
    f3.write("import os\n")
    f3.write("dir_list = [")
    f3.write(dir_list)
    f3.write("]\n")
    f3.write("for i in dir_list:\n")
    f3.write("    i = str(i)\n")
    f3.write("    os.chdir(i)\n")
    f3.write("    os.system('sbatch job.sh')\n")
'''        
    
        
        
    

