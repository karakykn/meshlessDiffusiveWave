import numpy as np
import matplotlib.pyplot as plt

def exact(x):
    k = 2
    return np.exp(k * x), k * np.exp(k *x), k**2 * np.exp(k *x)


def genMQ(x, sP):
    n = len(x)
    r = np.zeros((n, n))
    xdif = np.zeros((n, n))
    f = np.zeros((n, n + 3))
    fx = np.zeros((n, n + 3))
    fxx = np.zeros((n, n + 3))

    for i in range(n):
        for j in range(n):
            xdif[i, j] = x[i] - x[j]

    r[:, :] = np.abs(xdif[:, :])
    f[:n, :n] = np.sqrt(r[:, :] ** 2 + sP ** 2)
    f[:n, n] = 1
    f[:n, n + 1] = x[:]
    f[:n, n + 2] = x[:] ** 2
    fx[:n, :n] = xdif[:, :] / f[:n, :n]
    fx[:n, n + 1] = 1
    fx[:n, n + 2] = 2 * x[:]
    fxx[:n, :n] = 1 / f[:n, :n] - xdif[:, :] ** 2 / f[:n, :n] ** 3
    fxx[:n, n + 2] = 2

    S = np.zeros((n + 3, n + 3))
    S[0, :] = f[0, :]
    S[1:n - 1, :] = c * fx[1:-1, :] - D * fxx[1:-1, :]
    S[n - 1, :] = fx[-1, :]

    S[n, :n] = 1
    S[n + 1, :n] = x[:]
    S[n + 2, :n] = x[:] ** 2
    invS = np.linalg.pinv(S)
    sys = np.matmul(f, invS)
    return f,fx,fxx,sys

k = 2
c, D = 1, 1 / k
n = 5
x = np.linspace(0, 1, n)
sP = 4 * (x[1] - x[2])

f,fx,fxx,sys = genMQ(x, sP)
rhs = np.zeros(n+3)
rhs[0], rhs[n-1] = exact(x[0])[0], exact(x[-1])[1]

u = np.matmul(sys, rhs)
plt.plot(x, exact(x)[0])
plt.plot(x, u)
# plt.show()

# n_d = int(np.floor(n/2)+1)
# xL, xR = np.zeros(int(np.floor(n/2)+1)), np.zeros(int(np.floor(n/2)+1))
# xL[:], xR[:] = x[:int(np.floor(n/2))+1], x[int(np.floor(n/2)):]
#
# fL,fxL,fxxL,sysL = genMQ(xL, sP)
# fR,fxR,fxxR,sysR = genMQ(xR, sP)
# S_d = np.zeros((n+4, n+4))
#
# S_d[0, 0:n_d] = fL[0, 0:n_d]
# S_d[0, -3:] = fL[0, -3:]
# S_d[1:n_d-1, 0:n_d] = sysL[1:n_d-1, 0:n_d]
# S_d[1:n_d-1, -3:] = sysL[1:n_d-1, -3:]
#
# S_d[n_d+1:, 0:n_d]

n_d = int(np.floor(n/2)+1)
xL, xR = np.zeros(int(np.floor(n/2)+1)), np.zeros(int(np.floor(n/2)+1))
xL[:], xR[:] = x[:int(np.floor(n/2))+1], x[int(np.floor(n/2)):]

fL,fxL,fxxL,sysL = genMQ(xL, sP)
fR,fxR,fxxR,sysR = genMQ(xR, sP)
S_d = np.zeros((n+1, n+1))

S_d[0, 0:n_d] = fL[0, 0:n_d]
S_d[1:n_d-1, 0:n_d] = sysL[1:n_d-1, 0:n_d]
# S_d[-1, n_d:] =


# S_d[n_d+1:, 0:n_d]

pass
