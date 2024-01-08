from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import BFGS

# 创建一个简单的分子，比如水分子
h2o = molecule('H2O')

# 使用局域密度近似 (LDA) 的计算器
calc = EMT()

# 将计算器与 ASE 结构关联
h2o.set_calculator(calc)

# 创建优化器并运行结构优化
optimizer = BFGS(h2o)
optimizer.run(fmax=0.02)

# from ase.io import read, write
# from ase.optimize import BFGS
# from ase.calculators.lj import LennardJones
#
# # 从CIF文件中读取结构
# structure = read('C:\\Users\\sy\\Desktop\\SnH5CI3N2.cif')
#
# # 设置Lennard-Jones势计算器
# lj_calc = LennardJones()
#
# # 将计算器关联到结构
# structure.set_calculator(lj_calc)
#
# # 定义优化器，这里使用BFGS优化器
# optimizer = BFGS(structure, trajectory='optimized.traj')
#
# # 进行结构优化
# optimizer.run(fmax=0.05)  # 设置优化的力的阈值
#
# # 将优化后的结构写入新的CIF文件
# write('C:\\Users\\sy\\Desktop\\optimized.cif', structure, format='cif')