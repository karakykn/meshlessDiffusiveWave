# Meshless diffusive wave model for flood routing in channel networks

Code and results of

> Karakan, I., Jamal, J. F., Goharian, E., Chaudhry, M. H., and Imran, J. "Meshless solution of the diffusive wave
> equation for flood routing in channel networks." *Journal of Hydraulic Engineering* (manuscript HYENG-15162,
> under review).

Contact: Ismet Karakan, karakan@sc.edu

The model solves the diffusive wave equation for the discharge with Kansa's radial basis function collocation
method (multiquadric basis) and the theta (Crank-Nicolson) time scheme, on a single channel or a dendritic channel
network. Flow depths follow from the zero-inertia momentum equation; junctions are coupled through the volume
balance of the channels draining into them. The equations and the solution procedure are summarized at the top of
`source.py`.

## Requirements

Python 3 with NumPy for the model (`source.py`, `process.py`). The figure, table and sensitivity scripts also need
Matplotlib and pandas.

## Running a case

From this folder:

    python process.py syntheticChannel          # single channel (verification)
    python process.py akanYenNetwork            # backwater network of Akan and Yen (1981) (verification)
    python process.py lowerMississippi_paper    # Lower Mississippi River, 2011 (validation)
    python process.py 2008flood_stLouis         # Midwest flood, 2008 (validation)

Results are written to `<case>/segment<i>/run/<time>/Q` and `h` (discharge and depth at the nodes). The folders in
this repository hold the results of these runs.

## Case folder

| File | Content |
| --- | --- |
| `input` | Courant number; basis type (1: multiquadric); shape parameter coefficient c / r_av (used without a `shape` file); end time (s); print step (time steps between written results); outlet condition (0: normal depth, 1: depth from `boundary_h`) |
| `network` | one row `upstream downstream` for every pair of channels joined at a junction (`0 -1` for a single channel) |
| `shape` | c / r_av of every channel (segment order), selected by `select_shape.py` |
| `segment<i>/geo/nodes`, `slopes`, `mannings_n` | node positions (m), bed slope and Manning's n at the nodes |
| `segment<i>/geo/xsInfo`, `xSecs/xs<k>` | cross-section number of every node; cross-sections as rows (height above the bed, left bank x, right bank x) |
| `segment<i>/geo/boundary_Q`, `boundary_h` | inflow hydrograph [t (s), Q (m3/s)] of a channel without upstream channels; depth [t (s), h (m)] at the outlet. For channels at a junction the model writes the junction values here. |
| `segment<i>/geo/laterals`, `lateralDatas/` | point lateral flows (name, x, file); files hold [t (s), Q (m3/s)], negative for an outflow |
| `segment<i>/run/0/Q`, `h` | initial discharge and depth |

## Paper figures and tables

| Item | Script | Data |
| --- | --- | --- |
| Table 2 (shape parameter) | `select_shape.py` | case `shape` files, `sensitivity/shape_selection/` |
| Table 4 (error metrics) | `paper_scripts/table_metrics.py` | stored hydrographs, `<case>/data/` |
| Fig. 2 | `syntheticChannel/scripts/plot_single.py` | HEC-RAS and CNS series in `syntheticChannel/data/` |
| Fig. 3, Fig. S8 | `process_B1.py`, `paper_scripts/make_B1_full.py` | `sensitivity/B1_final_model/` |
| Fig. 4 | `process_B2.py` | `sensitivity/B2_final_model/` |
| Fig. 5 | `akanYenNetwork/scripts/hydrograph_figure.py` | `akanYenNetwork/data/` |
| Figs. 7 and 8 | `lowerMississippi_paper/scripts/boundary_figure.py`, `hydrograph_figure.py` | `lowerMississippi_paper/data/` |
| Figs. 10 and 11 | `2008flood_stLouis/scripts/boundary_figure.py`, `hydrograph_figure.py` | `2008flood_stLouis/data/` |
| Figs. S6 and S7 | `paper_scripts/fig_lateral_response.py`, `paper_scripts/lmr_profiles.py`, `fig_lateral_lmr.py` | `paper_scripts/` |
| CPU times of the LMR and Midwest cases (Table 4) | `time_meshless.py` | `timing_meshless.csv` |

The hydrograph scripts use the stored model series (`<case>/scripts/hydrographs_meshless.csv`) and rerun the model
with `--run`. The sweep scripts redraw their figures from the stored results with `--plot-only`; without it they
run the sweep (`--help` lists the options). `figure_tools.py` and `paper_style.py` hold the common figure code and
style. `akanYenNetwork/scripts/make_akanYen.py` builds the case folder of the backwater network. CNS results
(Beg et al. 2023, and the CNS runs on the inputs of these cases) are included as data; the CNS code is not part of
this repository, nor are the timing runs of Fig. 12 and the CPU times of the verification cases, which were run
together with CNS.
