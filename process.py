from source import *

caseName = 'syntheticNetwork_fine'
# caseName = 'syntheticChannel'
# caseName = 'lowerMississippi'
# caseName = '2008flood_stLouis'

solver = Network(caseName)
# solver.warmup(tol = 1e-8)
solver.solve()

'''For single channel, self.solve_ydk() is more accurate due to direct
interpolation of flow depths.'''
# solver.solve_ydk()
''''''
