"""Response of the discharge to a point lateral outflow (SM, Representation of Lateral Outflows).
Constant C and D, long channel, Q = Q0 before a lateral outflow Q_l < 0 at x0 starts at t = 0 (step):
    G = Q_l H(x - x0),   Q - G = Q0 - (Q_l / 2) erfc[(x0 + C t - x) / (2 sqrt(D t))]   (equation for Q - G, SM)
Steady state with C q_l alone (without -D dq_l/dx): Q = Q0 + Q_l exp[C (x - x0) / D] for x < x0, Q0 + Q_l for x > x0.
Dimensionless: xi = (x - x0) C / D, tau = C^2 t / D, discharges (. - Q0) / |Q_l|.
    python3 fig_lateral_response.py [--out folder]
    python3 fig_lateral_response.py [--out folder]                               (Mac, Helvetica)
Output: lateral_response.pdf/.png (default: this folder)
"""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_style import PAPER_STYLE, N_RAMP, REFERENCE_GREY, savefig_safe

out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else HERE
plt.rcdefaults()
plt.rcParams.update(PAPER_STYLE)
erfc = np.vectorize(math.erfc)

xi = np.linspace(-4, 12, 3201)
up, dn = xi[xi < 0], xi[xi > 0]
TAUS = [0.0, 1.0, 4.0, np.inf]
LAB = [r'$C^2t/D=0^+$', r'$C^2t/D=1$', r'$C^2t/D=4$', r'$C^2t/D\rightarrow\infty$']
STY = [dict(color=N_RAMP[0], ls=':'), dict(color=N_RAMP[1], ls='-.'), dict(color=N_RAMP[2], ls='--'),
       dict(color=N_RAMP[3], ls='-')]


def P(x, tau):
    """(P - Q0)/|Q_l|"""
    if tau == 0:
        return (x > 0).astype(float)
    if np.isinf(tau):
        return np.zeros_like(x)
    return 0.5 * erfc((tau - x) / (2 * math.sqrt(tau)))


fig, (a, b) = plt.subplots(2, 1, figsize=(3.5, 5.1), sharex=True)
for ax in (a, b):
    ax.axvline(0, color=REFERENCE_GREY, lw=0.8, zorder=0)
handles = []
# (a) G and P
a.plot(dn, -np.ones_like(dn), color='black', lw=1.0)
a.plot([xi[0], 0], [0, 0], color='black', lw=1.0)
a.text(6.0, -0.9, r'$G$, all $t>0$', ha='center', va='bottom', fontsize=8)
for tau, lab, st in zip(TAUS, LAB, STY):
    if tau == 0:
        h, = a.plot(dn, P(dn, tau), lw=1.5, **st, label=lab)
        a.plot(up, P(up, tau), lw=1.5, **st)
    else:
        h, = a.plot(xi, P(xi, tau), lw=1.5, **st, label=lab)
    handles.append(h)
a.text(6.0, 1.06, r'$Q-G$', ha='center', va='bottom', fontsize=8)
a.set_ylim(-1.25, 1.35)
a.set_ylabel(r'$(Q-G-Q_0)/|Q_\ell|$,  $G/|Q_\ell|$')
a.text(0.15, -1.17, r'$x_0$', ha='left', va='bottom', fontsize=8)
# (b) Q = P + G, and the steady state with C q_l alone
g, = b.plot(up, -np.exp(up), color=REFERENCE_GREY, lw=3.0, solid_capstyle='butt', label=r'$C^2t/D\rightarrow\infty$, $Cq_\ell$ alone')
b.plot(dn, -np.ones_like(dn), color=REFERENCE_GREY, lw=3.0, solid_capstyle='butt')
for tau, lab, st in zip(TAUS, LAB, STY):
    b.plot(up, P(up, tau), lw=1.5, **st)
    b.plot(dn, P(dn, tau) - 1.0, lw=1.5, **st)
handles.append(g)
b.set_ylim(-1.25, 0.5)
b.set_xlim(xi[0], xi[-1])
b.set_xticks(np.arange(-4, 13, 2))
b.set_ylabel(r'$(Q-Q_0)/|Q_\ell|$')
b.set_xlabel(r'$(x-x_0)\,C/D$')
for ax, l in ((a, 'a'), (b, 'b')):
    ax.text(0.03, 0.95, f'({l})', transform=ax.transAxes, fontweight='bold', va='top', ha='left')
fig.tight_layout(h_pad=0.6, rect=(0, 0, 1, 0.885))
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.56, 1.0), ncol=2, fontsize=8,
           handlelength=2.6, columnspacing=1.0)
os.makedirs(out, exist_ok=True)
for ext in ('pdf', 'png'):
    savefig_safe(fig, os.path.join(out, 'lateral_response.' + ext), bbox_inches='tight')
print('Saved', os.path.join(out, 'lateral_response.pdf'))
