"""Runs a case with the meshless model (source.py) and writes the results to <case>/segment<i>/run/<time>/Q and h,
every 'print step' time steps of the input file and at the end time.

    python process.py syntheticChannel          single channel (verification)
    python process.py akanYenNetwork            backwater network of Akan and Yen (1981) (verification)
    python process.py lowerMississippi_paper    Lower Mississippi River, 2011 (validation)
    python process.py 2008flood_stLouis         Midwest flood, 2008 (validation)
"""
import sys

from source import Network

THETA = 0.5     # weight of the new time level (0.5: Crank-Nicolson), all cases of the paper

if __name__ == '__main__':
    case = sys.argv[1] if len(sys.argv) > 1 else 'syntheticChannel'
    Network(case).solve(THETA)
