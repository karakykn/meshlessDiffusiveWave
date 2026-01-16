import numpy as np
import matplotlib.pyplot as plt

def genMQ(x, c=0):
    n = len(x)
    xdif = np.zeros((n,n))
    r = np.zeros((n, n))
    f = np.zeros((n,n))
    fx = np.zeros((n,n))
    fxx = np.zeros((n,n))
    rav = 0
    for i in range(n):
        for j in range(n):
            xdif[i,j] = x[i] - x[j]
    for i in range(n-1):
        rav += (x[i+1] - x[i]) / n

    c = 4 * rav
    r[:,:] = np.abs(xdif[:,:])
    f[:,:] = np.sqrt(r[:,:]**2 + c**2)
    fx[:,:] = xdif[:,:] / f[:,:]
    fxx[:,:] = 1 / f[:, :] - xdif[:, :] ** 2 / f[:, :] ** 3
    return f,fx,fxx

def solve_diffusion(x, dBc, nBc, split, c=0):
    f,fx,fxx = genMQ(x, c)
    n = fxx.shape[0]
    sys = np.zeros((n,n))
    rhs = np.zeros(n)
    sys[:,:] = fxx[:,:]
    sys[0,:] = f[0,:]
    sys[-1,:] = fx[-1,:]
    rhs[0] = dBc
    rhs[-1] = nBc
    sys, rhs, f = decompose_sys_diff(sys, rhs, split, fx, f)
    invSys = np.linalg.pinv(sys)
    ind = np.linspace(0,n, n+1)
    mask = ind != split
    soln = np.matmul(f, np.matmul(invSys, rhs))[mask]
    return soln

def decompose_sys_diff(sys, rhs, split, fx, f):
    n = sys.shape[0]
    sys_n = np.zeros((n + 1, n + 1))
    rhs_n = np.zeros(n + 1)
    f_n = np.zeros((n + 1, n + 1))
    sys_n[:split,:split+1] = sys[:split,:split+1]
    rhs_n[:split] = rhs[:split]

    rhs_n[split+2:] = rhs[split+1:]
    sys_n[split+2:, split+1:] = sys[split+1:, split:]

    sys_n[split, :split+1] = fx[split,:split+1]
    sys_n[split, split+1:] = -fx[split,split:]

    sys_n[split + 1, :split+1] = f[split,:split+1]
    sys_n[split + 1, split+1:] = -f[split,split:]

    f_n[:split+1, :split+1] = f[:split+1, :split+1]
    f_n[split + 1:, split + 1:] = f[split:, split:]

    return sys_n, rhs_n, f_n

def solve_diffusion_fdm(x, dBc, nBc, split):
    n = len(x)
    sys = np.zeros((n,n))
    rhs = np.zeros(n)
    dx = (x[-1] - x[0]) / (n-1)
    for i in range(1, n-1):
        sys[i, i] = -2 / dx**2
        sys[i, i+1] = 1 / dx**2
        sys[i, i-1] = 1 / dx**2
    sys[0,0] = 1
    sys[-1,-1] = 1 / dx
    sys[-1, -2] = -1 / dx
    rhs[0] = dBc
    rhs[-1] = nBc
    sys, rhs = decompose_sys_diffusion(sys, rhs, split)
    invSys = np.linalg.pinv(sys)
    ind = np.linspace(0,n, n+1)
    mask = ind != split
    soln = np.matmul(invSys, rhs)[mask]
    return soln

def decompose_sys_diffusion(sys, rhs, split):
    n = sys.shape[0]
    sys_n = np.zeros((n + 1, n + 1))
    rhs_n = np.zeros(n + 1)
    sys_n[:split,:split+1] = sys[:split,:split+1]
    rhs_n[:split] = rhs[:split]

    rhs_n[split+2:] = rhs[split+1:]
    sys_n[split+2:, split+1:] = sys[split+1:, split:]

    sys_n[split + 1, split + 1] = -1
    sys_n[split + 1, split] = 1

    sys_n[split, split -1:split+1] = sys[-1, -2:]
    sys_n[split, split+1: split+3] = -sys[-1, -2:]

    return sys_n, rhs_n

n = 5
split = int((n-1)/2)
x = np.linspace(0, 1, n)

solution = solve_diffusion(x, 1, 1, split)
# solution = solve_diffusion_fdm(x, 1, 1, split)
fig, ax1 = plt.subplots(1, figsize=(7,4))
ax1.plot(x, solution)
plt.show()