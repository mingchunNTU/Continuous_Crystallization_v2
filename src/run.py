from data import *
from crystal import *

# ===========================================

setting_dir=r"../examples/ammonium alum/"
CSD_mesh_size=5 # um

# ===========================================

file=setting_dir+"setting.csv"
tmp1=crystal()
tmp1.read_setting(file)

tmp1.solve()
tmp1.output(setting_dir,CSD_mesh_size)



