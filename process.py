from source import *

# caseName = 'syntheticNetwork'
# caseName = 'syntheticChannel'
# caseName = 'lowerMississippi'
# caseName = '2008flood_stLouis'
caseName = 'folsom/2017'
# caseName = 'folsom/2017_single'

solver = Network(caseName)
# solver.warmup(tol = 1e-3)
solver.solve()

# solver.warmup2(tol = 1e-3)
# solver.solve2()


