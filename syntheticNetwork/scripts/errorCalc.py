import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def rmse(exact, approx):
    return np.sqrt(np.sum((exact-approx) ** 2) / len(exact))

def mpe(exact, approx):
    return np.sum(np.abs((exact - approx) / exact)) / len(exact) * 100

caseName = '../'
hecras_dir = os.path.join('..', 'data', 'ex2')
nwm_dir = os.path.join('..', 'data', 'ex2_nwm')