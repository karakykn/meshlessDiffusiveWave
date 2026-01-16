from source import *

# caseName = 'syntheticNetwork'
# caseName = 'syntheticChannel'
# caseName = 'lowerMississippi'
caseName = '2008flood_stLouis'

solver = Network(caseName)
# solver.warmup(tol = 1e-3)
solver.solve()


