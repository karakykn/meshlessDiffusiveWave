"""Figure style shared by every figure of the paper: the hydrographs (figure_tools.py and
<case>/scripts/hydrograph_figure.py) and the sensitivity sweeps (process_B1.py, process_B2.py).

Helvetica, boxed axes with inward ticks, no grid, a boxed legend. Series are told apart by a blue that gets
darker (lightness L* about 72, 55, 39 and 23, so they print as evenly spaced greys), a line style and a
marker; references (HEC-RAS, the parameter values of the cases, zero lines) are light grey."""

PAPER_STYLE = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'sans', 'mathtext.it': 'sans:italic', 'mathtext.default': 'regular',
    'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8.5,
    'axes.linewidth': 0.8, 'axes.edgecolor': 'black', 'axes.labelcolor': 'black', 'text.color': 'black',
    'xtick.color': 'black', 'ytick.color': 'black',
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 3.5, 'ytick.major.size': 3.5, 'xtick.minor.size': 2, 'ytick.minor.size': 2,
    'lines.linewidth': 1.3, 'lines.markersize': 3.5, 'lines.markeredgewidth': 0.8,
    'legend.frameon': True, 'legend.edgecolor': 'black', 'legend.fancybox': False, 'legend.framealpha': 1,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300,
}

# Ordered series (N, Courant number, ...): light to dark blue, each with its own line style and marker.
N_RAMP = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
N_STYLES = [dict(ls=':', marker='o', filled=False), dict(ls='-.', marker='^', filled=False),
            dict(ls='--', marker='s', filled=True), dict(ls='-', marker='o', filled=True),
            dict(ls=(0, (6, 2, 1, 2, 1, 2)), marker='D', filled=False), dict(ls='-', marker='v', filled=False)]
REFERENCE_GREY = '0.7'                     # reference lines and the HEC-RAS band

# Hydrographs: the meshless model takes the darkest blue (solid, filled circles) and CNS the second blue
# (dash-dot, open triangles), as in the sweeps; HEC-RAS is a wide grey band under them.
HYDRO_SERIES = {
    'mesh':   dict(color=N_RAMP[3], mec=N_RAMP[3], mfc=N_RAMP[3], ls='-', lw=1.5, marker='o', label='Meshless', zorder=5),
    'cns':    dict(color=N_RAMP[1], mec=N_RAMP[1], mfc='white', ls='-.', lw=1.3, marker='^', ms=4.0,
                   label='CNS (Beg et al. 2023)', zorder=4),
    'hec':    dict(color=REFERENCE_GREY, ls='-', lw=3.4, solid_capstyle='butt', label='Dynamic (HEC-RAS)', zorder=2),
    'inflow': dict(color='black', ls=':', lw=0.9, label='Upstream inflow', zorder=1),
    'obs':    dict(color='black', ls='none', marker='o', ms=2.4, mfc='none', mew=0.5, label='USGS', zorder=3),
}


def savefig_safe(fig, path, **kw):
    """fig.savefig(path, **kw). If a PDF cannot be written with the fonts embedded as TrueType (Type 42), which
    some matplotlib versions fail to do with the macOS .ttc fonts (Helvetica), the error is printed and the PDF is
    written with Type 3 fonts instead."""
    import matplotlib.pyplot as plt
    try:
        fig.savefig(path, **kw)
    except Exception as e:
        if not str(path).lower().endswith('.pdf'):
            raise
        print(f'Saving {path} with TrueType fonts failed ({type(e).__name__}: {e}); writing it with Type 3 fonts.')
        with plt.rc_context({'pdf.fonttype': 3}):
            fig.savefig(path, **kw)
