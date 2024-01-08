from ase.io import write
from ase.build import molecule
from ase.calculators.lammpslib import LAMMPSlib
from ase.md import VelocityVerlet

# 创建一个简单的分子，比如水分子
atoms = molecule('H2O')

# 将 LAMMPS 设置为 ASE 的计算器
lammps_cmd = 'lmp_serial'  # 替换为你的 LAMMPS 可执行文件路径
parameters = {'pair_style': 'lj/cut/coul/long 10.0',
              'pair_coeff': ['1 1 1.0 1.0', '2 2 1.0 1.0', '3 3 1.0 1.0'],
              'mass': ['1 1.0', '2 1.0', '3 1.0'],
              'boundary': 'p p p'}
calc = LAMMPSlib(lmpcmds=parameters, lammpsrun=lammps_cmd, log_file='lammps.log')

# 将计算器与 ASE 结构关联
atoms.set_calculator(calc)

# 进行分子动力学模拟
dyn = VelocityVerlet(atoms, timestep=1.0)  # 时间步长可根据需要调整
dyn.run(100)  # 运行100个步长的分子动力学模拟
