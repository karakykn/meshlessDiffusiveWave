"""
Meshless (RBF collocation) diffusive wave model for a single channel or a dendritic channel network.

Discharge: Q-form of the diffusive wave equation with lateral flow (Cappelaere 1997),
        dQ/dt + C dQ/dx = D d2Q/dx2 + C q_l - D dq_l/dx.
A lateral flow Q_l,k(t) at x_k (segment<i>/geo/laterals and lateralDatas; negative for an outflow such as a
spillway) is a point source, q_l = sum_k Q_l,k delta(x - x_k). With
        Q = P + G,   G(x, t) = sum_k Q_l,k(t) W_k(x),   W_k = 0 upstream of x_k, 1 downstream (1/2 at a node on x_k),
the equation becomes dP/dt + C dP/dx = D d2P/dx2 - dG/dt. P is solved with the collocation system and Q = P + G,
so the lateral flow is never differentiated and leaves the channel exactly.

Parameters: Qbar, Cbar and Dbar are the uniform-flow values at the local depth h (tabulated per node) and
        Gamma = Q / Qbar(h)
        C = Cbar / 2 [Gamma (1 + k) + (1 - k) / Gamma],   k = (Qbar / Dbar) dDbar/dQbar   (Cappelaere 1997, Eq. 14a)
        D = Dbar / |Gamma|.
Depth: zero-inertia momentum equation, integrated upstream from the downstream depth,
        dh/dx = S0 - S_f = S0 (1 - Gamma |Gamma|).

Spatial discretization: Kansa's collocation method with the multiquadric (MQ) basis, sqrt(r^2 + c^2), plus a
quadratic polynomial; c = coef * r_av, r_av = L / (N - 1) the mean node spacing. coef is given per channel in the
file <case>/shape (one value per segment, written by select_shape_coefficients) or, without that file, by the
input file for all channels. Time discretization: theta scheme (theta = 0.5: Crank-Nicolson).

One time step (Network.step):
 1. C and D at every node from the current Q and h; dt = Courant number * min(dx / C) over the network.
 2. Forward (discharge) sequence, from the upstream channels down: the inflow at the upstream end (hydrograph, or
    the sum of the outflows of the channels draining into the junction) and, at the downstream end, the outflow
    row dQ/dt + C dQ/dx = 0 at the network outlet, or the junction outflow Q_J of a channel draining into a
    junction.
 3. Junction outflow (networks only): Q_J of every channel draining into a junction satisfies the volume balance
    of the channel,
        S^{n+1} - S^n = dt [theta (Q_in - Q_J)^{n+1} + (1 - theta) (Q_in - Q_J)^n] + lateral volume,
    S = integral of the wetted area A(h) dx. Within a step Q is linear in Q_in and Q_J, and S and the depth at the
    upstream end of every channel are linearized with the Jacobian of the depth system of the previous step, so all
    balances form one small linear system for the Q_J. The residual of the linearization is carried into the
    balance of the next step.
 4. Backward (depth) sequence, from the outlet up: the depth at the downstream end (rating curve or prescribed
    depth at the outlet; at a junction, the depth at the upstream end of the downstream channel), the depths along
    the channel from the momentum equation (Newton's method), then Gamma = Q / Qbar(h).

Case folder (see README.md): input, network, optional shape, and segment<i>/geo (nodes, slopes, mannings_n, xsInfo,
xSecs/xs<k>, boundary_Q, boundary_h, optional laterals and lateralDatas) and segment<i>/run/0 (initial Q and h).
Results are written to segment<i>/run/<time>/Q and h.
"""
import csv
import math
import os
import shutil
from collections import defaultdict, deque
from time import perf_counter, process_time

import numpy as np

try:            # LAPACK gesv without the argument checks of np.linalg.solve (same routine, identical results)
    from numpy.linalg import _umath_linalg as _UL

    def _solve(a, b):
        return _UL.solve(a, b, signature='dd->d')
    _solve(np.eye(2), np.ones((2, 1)))
except Exception:
    _solve = np.linalg.solve

N_TABLE = 300           # intervals of the depth grid of the cross-section and uniform-flow tables


# ---------------------------------------------------------------------------------------------------------------------
# Table lookups
# ---------------------------------------------------------------------------------------------------------------------
def _interp_table(x, tab, rows):
    """Row-wise, vectorized np.interp: for every node i, np.interp(x[i], tab[i, 0, :], tab[i, r, :]) for each r in
    rows, as an array of shape (len(rows), n). Same arithmetic as np.interp (clamping at both ends, exact hits on
    grid points, NaN handling)."""
    xp = tab[:, 0, :]
    n, m = xp.shape
    node = np.arange(n)
    j = np.count_nonzero(xp <= x[:, None], axis=1) - 1      # last grid point <= x (-1: below the table)
    jc = np.clip(j, 0, m - 2)
    x0 = xp[node, jc]
    x1 = xp[node, jc + 1]
    r = np.asarray(rows)[:, None]
    y0 = tab[node, r, jc]                                     # shape (len(rows), n)
    y1 = tab[node, r, jc + 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = (y1 - y0) / (x1 - x0)
        out = slope * (x - x0) + y0
        bad = np.isnan(out)
        if bad.any():                                         # same fallback as numpy
            out = np.where(bad, slope * (x - x1) + y1, out)
            out = np.where(bad & np.isnan(out) & (y0 == y1), y0, out)
    out = np.where(x0 == x, y0, out)                          # exactly on a grid point
    out = np.where(j < 0, tab[node, r, 0], out)               # below the table
    out = np.where(j >= m - 1, tab[node, r, m - 1], out)      # at / above the top of the table
    return np.where(np.isnan(x), x, out)


class _NodeTable(object):
    """Lookup tables of all nodes of a channel on one uniform depth grid per node (rows: Qbar, Cbar, Dbar,
    dDbar/dQbar, top width, area). values() returns the same numbers as _interp_table, but finds the grid interval
    from h / dh instead of comparing h with the whole grid, and one interval search serves all rows."""

    def __init__(self, grid, values):
        self.g = np.ascontiguousarray(grid, dtype=float)                      # (N, m)
        self.v = np.ascontiguousarray(values, dtype=float)                    # (N, R, m)
        n, m = self.g.shape
        self.m = m
        self.node = np.arange(n)
        self.shared = bool(np.all(self.g == self.g[0]))                      # one grid for all nodes
        self.g1 = self.g[0].copy()
        with np.errstate(divide='ignore', invalid='ignore'):
            self.inv_dh = (m - 1) / (self.g[:, -1] - self.g[:, 0])
            self.slope = (self.v[:, :, 1:] - self.v[:, :, :-1]) / (self.g[:, None, 1:] - self.g[:, None, :-1])
            self.ok = bool(np.all(self.g[:, 0] == 0) and np.all(np.isfinite(self.inv_dh)) and
                           np.allclose(np.diff(self.g, axis=1) * self.inv_dh[:, None], 1, rtol=1e-9, atol=0))
        if self.shared:
            self.inv_dh = float(self.inv_dh[0])
        # rows whose values and slopes are finite everywhere need no NaN fallback
        self.finite = np.all(np.isfinite(self.v), axis=(0, 2)) & np.all(np.isfinite(self.slope), axis=(0, 2))
        self._rows = {}
        if self.ok:                         # the interval search must agree with _interp_table exactly
            xs = np.concatenate([self.g.T, (self.g[:, 1:] + self.g[:, :-1]).T / 2, self.g.T * (1 + 1e-15),
                                 self.g.T * (1 - 1e-15), -self.g[:, 1:2].T, 2 * self.g[:, -1:].T])
            for x in xs:
                if not np.array_equal(self._interval(x), np.count_nonzero(self.g <= x[:, None], axis=1) - 1):
                    self.ok = False
                    break

    def _interval(self, x):
        """index of the last grid point <= x (-1 below the table)"""
        jc = (x * self.inv_dh).astype(np.intp)            # = floor for x >= 0 (index() passes finite x only)
        np.minimum(jc, self.m - 2, out=jc)
        np.maximum(jc, 0, out=jc)
        if self.shared:
            return jc - (self.g1[jc] > x) + (self.g1[jc + 1] <= x)
        return jc - (self.g[self.node, jc] > x) + (self.g[self.node, jc + 1] <= x)

    def index(self, x):
        """interval search for the depths x (None for a NaN or infinite depth: _interp_table handles it)"""
        if not np.isfinite(x).all():
            return None
        j = self._interval(x)
        jmin, jmax = j.min(), j.max()
        jc = j if (jmin >= 0 and jmax <= self.m - 2) else np.clip(j, 0, self.m - 2)
        x0 = self.g1[jc] if self.shared else self.g[self.node, jc]
        return x, j, jc, x0, jmin, jmax

    def values(self, state, rows, x=None):
        """table values of the given rows at the depths of an index() result"""
        if state is None:
            tab = np.concatenate([self.g[:, None, :], self.v], axis=1)
            return _interp_table(x, tab, np.asarray(rows) + 1)
        x, j, jc, x0, jmin, jmax = state
        r = self._rows.get(rows)
        if r is None:
            r = self._rows[rows] = (np.asarray(rows)[:, None], bool(self.finite[list(rows)].all()))
        r, finite = r
        y0 = self.v[self.node, r, jc]
        slope = self.slope[self.node, r, jc]
        out = slope * (x - x0) + y0
        if not finite:
            bad = np.isnan(out)
            if bad.any():                                     # same fallback as numpy
                x1 = self.g[self.node, jc + 1]
                y1 = self.v[self.node, r, jc + 1]
                with np.errstate(invalid='ignore'):
                    out = np.where(bad, slope * (x - x1) + y1, out)
                out = np.where(bad & np.isnan(out) & (y0 == y1), y0, out)
            hit = x0 == x
            if hit.any():                                     # exactly on a grid point
                out = np.where(hit, y0, out)
        if jmin < 0:                                          # below the table
            out = np.where(j < 0, self.v[self.node, r, 0], out)
        if jmax >= self.m - 1:                                # at / above the top of the table
            out = np.where(j >= self.m - 1, self.v[self.node, r, self.m - 1], out)
        return out


def _write_column(path, arr):
    """same bytes as np.savetxt(path, arr) for a 1-D float array"""
    with open(path, 'w') as fh:
        fh.write(''.join(['%.18e\n' % v for v in arr]))


class _JunctionSeries(object):
    """In-memory record [time, value] of a junction boundary: the inflow of a channel that starts at a junction (sum
    of the outflows of the channels draining into it) or the downstream depth of a channel that ends at a junction.
    Network mirrors the records to segment<j>/geo/boundary_Q and boundary_h at every print step."""
    __slots__ = ('t', 'v', '_pos')

    def __init__(self):
        self.t, self.v, self._pos = [], [], {}

    def _append(self, time, val):
        self._pos[time] = len(self.t)
        self.t.append(time)
        self.v.append(val)

    def add(self, time, val):
        """value(time) += val, or a new row if the time is not recorded yet"""
        k = self._pos.get(time)
        if k is None:
            self._append(time, val)
        else:
            self.v[k] = self.v[k] + val

    def set(self, time, val):
        """value(time) = val, or a new row if the time is not recorded yet"""
        k = self._pos.get(time)
        if k is None:
            self._append(time, val)
        else:
            self.v[k] = val

    def value_at(self, time):
        """np.interp over the record; an exactly recorded time is returned directly"""
        k = self._pos.get(time)
        if k is not None:
            return self.v[k]
        return np.interp(time, np.asarray(self.t, dtype=float), np.asarray(self.v, dtype=float))

    def trim(self, n):
        """drop the last n rows but keep the latest one (keeps the record short during long runs)"""
        last_t, last_v = self.t[-1], self.v[-1]
        self.t = self.t[:-n] + [last_t]
        self.v = self.v[:-n] + [last_v]
        self._pos = {t: k for k, t in enumerate(self.t)}

    def as_array(self):
        return np.column_stack([np.asarray(self.t, dtype=float), np.asarray(self.v, dtype=float)])


# ---------------------------------------------------------------------------------------------------------------------
# Shape parameter
# ---------------------------------------------------------------------------------------------------------------------
def select_shape_coefficients(caseName, candidates=(2, 3, 4, 6, 8, 12, 16, 24, 32, 48), write=True):
    """Shape parameter coefficient (c = coef * r_av) of every channel of a case, chosen on the channel's own nodes:
    the candidate with the smallest error of the operators of the scheme, dQ/dx (fx_invF), d2Q/dx2 (fxx_invF, at the
    interior nodes) and the depth integration of the backward sequence (hsys), on two smooth test profiles over the
    channel, exp(-3(1 - s)) and sin(2 pi s) + s^2, s = (x - x_0) / L (largest relative errors). Writes
    <caseName>/shape (read by Network) and returns (coefficients, errors)."""
    import types
    segs = sorted((d for d in os.listdir(caseName) if d.startswith('segment')), key=lambda d: int(d[7:]))
    coefs, table = [], []
    for sname in segs:
        x = np.loadtxt(os.path.join(caseName, sname, 'geo', 'nodes'))
        N, L = len(x), x[-1] - x[0]
        xi = (x - x[0]) / L
        tests = ((np.exp(-3 * (1 - xi)), 3 * np.exp(-3 * (1 - xi)) / L, 9 * np.exp(-3 * (1 - xi)) / L ** 2),
                 (np.sin(2 * np.pi * xi) + xi ** 2, (2 * np.pi * np.cos(2 * np.pi * xi) + 2 * xi) / L,
                  (2 - 4 * np.pi ** 2 * np.sin(2 * np.pi * xi)) / L ** 2))
        errs = []
        for b in candidates:
            o = types.SimpleNamespace(nodeNo=N, geo={'nodes': x})
            o.f = np.zeros((N, N + 3))
            o.fx = np.zeros_like(o.f)
            o.fxx = np.zeros_like(o.f)
            SingleChannel.buildMQ(o, b)
            e = 0.0
            for F, dF, d2F in tests:
                e = max(e, np.max(np.abs(o.fx_invF @ F - dF)) / np.max(np.abs(dF)))
                e = max(e, np.max(np.abs((o.fxx_invF @ F - d2F)[1:-1])) / np.max(np.abs(d2F[1:-1])))
                R = np.zeros(N + 3)
                R[:N - 1] = dF[:-1]
                R[N - 1] = F[-1]
                e = max(e, np.max(np.abs(o.hsys @ R - F)) / np.max(np.abs(F)))
            errs.append(e)
        coefs.append(float(candidates[int(np.argmin(errs))]))
        table.append(errs)
    if write:
        np.savetxt(os.path.join(caseName, 'shape'), np.array(coefs), fmt='%g')
    return coefs, table


# ---------------------------------------------------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------------------------------------------------
class Network(object):
    """A case: one channel, or a dendritic network of channels (segments) joined at junctions.

    Input file (<case>/input, one value under each of six headers): Courant number, basis type (1: MQ, the only
    one implemented), shape parameter coefficient (used when <case>/shape does not exist), end time (s), print
    step (number of time steps between written results) and downstream boundary condition at the outlet
    (0: normal depth from the rating curve, 1: depth from boundary_h).
    Network file (<case>/network): one row 'upstream downstream' per pair of segments joined at a junction."""

    def __init__(self, caseName):
        self.caseName = caseName
        self.clean_run_directories()
        self.network = np.atleast_2d(np.loadtxt(caseName + '/network', dtype=int))   # one row: two channels

        with open(caseName + '/input', 'r') as file:
            lines = file.readlines()
        self.CFL = float(lines[1].strip())               # Courant number
        rbf_type = int(lines[4].strip())
        if rbf_type != 1:
            raise ValueError(f'{caseName}/input: basis type {rbf_type}; only the MQ basis (1) is implemented')
        self.shapeCoef = float(lines[7].strip())
        self.simEndTime = float(lines[10].strip())
        self.printStep = int(lines[13].strip())
        self.dsbcSpec = int(lines[16].strip())           # outlet: 0 normal depth, 1 depth from boundary_h
        self.time = 0

        self.num_segments = sum(1 for d in os.listdir(caseName)
                                if d.startswith('segment') and os.path.isdir(os.path.join(caseName, d)))
        self.segments = self.load_segments()

        # calcOrder: segments from upstream to downstream; segUpstreamInfo[i] / segDownstreamInfo[i]: the segments
        # that drain into segment i / that segment i drains into
        if self.num_segments == 1:
            self.calcOrder = [0]
            self.segUpstreamInfo = defaultdict(list)
            self.segDownstreamInfo = defaultdict(list)
        else:
            self.order_segments()

        self.reset_junctions()
        self._has_junctions = any(len(self.segDownstreamInfo[i]) > 0 for i in self.calcOrder)
        if self._has_junctions:
            for s in self.segments.values():
                s.prepare_linear()

    def clean_run_directories(self):
        """removes the results of a previous run: all time folders in segment<i>/run except '0'"""
        for segment_dir in os.listdir(self.caseName):
            run_dir = os.path.join(self.caseName, segment_dir, 'run')
            if segment_dir.startswith('segment') and os.path.isdir(run_dir):
                for time_dir in os.listdir(run_dir):
                    time_dir_path = os.path.join(run_dir, time_dir)
                    if os.path.isdir(time_dir_path) and time_dir != '0':
                        shutil.rmtree(time_dir_path)
        print(f"Completed cleaning up the 'run' directories in {self.caseName}.")

    def load_segments(self):
        """SingleChannel objects of all segments, with c = coef * r_av; coef from <case>/shape (one value per
        segment, segment order) if the file exists, otherwise the value of the input file for all segments"""
        coefs = [self.shapeCoef] * self.num_segments
        shape_file = os.path.join(self.caseName, 'shape')
        if os.path.exists(shape_file):
            v = np.atleast_1d(np.loadtxt(shape_file, dtype=float))
            if len(v) != self.num_segments:
                raise ValueError(f'{shape_file}: {len(v)} values for {self.num_segments} segments')
            coefs = [float(x) for x in v]
        self.shapeCoefs = coefs
        return {i: SingleChannel(self.caseName, i, coefs[i]) for i in range(self.num_segments)}

    def order_segments(self):
        """calcOrder (topological order, upstream segments first) and the junction connections from the network file"""
        self.segUpstreamInfo = defaultdict(list)
        self.segDownstreamInfo = defaultdict(list)
        in_degree = defaultdict(int)                     # number of segments draining into each segment
        all_nodes = set()
        for upstream, downstream in self.network:
            self.segUpstreamInfo[downstream].append(upstream)
            self.segDownstreamInfo[upstream].append(downstream)
            in_degree[downstream] += 1
            all_nodes.update([upstream, downstream])
        queue = deque([node for node in all_nodes if in_degree[node] == 0])
        order = []
        while queue:
            segment = queue.popleft()
            order.append(segment)
            for downstream in self.segDownstreamInfo[segment]:
                in_degree[downstream] -= 1
                if in_degree[downstream] == 0:
                    queue.append(downstream)
        self.calcOrder = np.array(order)

    # ---- time integration ------------------------------------------------------------------------------------------
    def solve(self, theta=0.5):
        """Runs the case to the end time of the input file (theta: weight of the new time level, 0.5 = Crank-Nicolson).
        Q and h of every channel are written to segment<i>/run/<time> every printStep steps and at the end."""
        nsteps = 0
        if self.time == 0:
            self.initial_depths()
        # CPU time of the time loop (process_time: CPU time of this process; perf_counter: wall clock); cpu_out: the
        # part spent on screen and file output every printStep steps
        cpu_out = 0.0
        cpu0, wall0 = process_time(), perf_counter()
        while self.time < self.simEndTime:
            self.step(theta)
            nsteps += 1
            if nsteps % self.printStep == 0:
                out0 = process_time()
                print(f'Time: {self.time:.2f} s')
                print('            max Q      min Q      max h      min h      max dQ/dx      min dQ/dx')
                for i in self.calcOrder:
                    s = self.segments[i]
                    print(f'Channel {i + 1}: {max(s.Q):10.2f} {min(s.Q):10.2f} {max(s.h):10.2f} {min(s.h):10.2f} '
                          f'{max(s.dQdx):14.6f} {min(s.dQdx):14.6f}')
                    self.write_state(i)
                self.trim_junction_records(self.printStep)
                print()
                cpu_out += process_time() - out0
        cpu_loop, wall_loop = process_time() - cpu0, perf_counter() - wall0
        for i in self.calcOrder:
            self.write_state(i)
        self._flush_junction_files()
        self.report_cpu_time(nsteps, cpu_loop, cpu_out, wall_loop)

    def report_cpu_time(self, nsteps, cpu_loop, cpu_out, wall_loop):
        """prints the CPU time of the time loop (also kept in self.cpu_time)"""
        nodes = sum(len(self.segments[i].Q) for i in range(self.num_segments))
        self.cpu_time = dict(case=self.caseName, channels=self.num_segments, nodes=nodes, steps=nsteps,
                             cpu_loop=cpu_loop, cpu_output=cpu_out, cpu_compute=cpu_loop - cpu_out,
                             wall_loop=wall_loop)
        print('================ CPU time of the time loop ================')
        print(f'case: {self.caseName} | channels: {self.num_segments} | nodes: {nodes} | time steps: {nsteps}')
        print(f'CPU time, time loop          : {cpu_loop:10.3f} s  ({1e3 * cpu_loop / max(nsteps, 1):.4f} ms/step)')
        print(f'  of which output (print/write): {cpu_out:8.3f} s')
        print(f'CPU time, without output     : {cpu_loop - cpu_out:10.3f} s')
        print(f'wall-clock time, time loop   : {wall_loop:10.3f} s')
        print('============================================================')

    def initial_depths(self):
        """replaces the initial depths (run/0/h) by the depths of the initial discharges (momentum equation from the
        downstream depth, as in the backward sequence); the discharges are not changed"""
        self._backward_sequence()
        if self._has_junctions:
            self._prepare_volume_junction()

    def _time_step(self):
        """C and D at every node from the current Q and h, and the time step of the whole network"""
        dt = 1e8
        for i in range(self.num_segments):
            seg = self.segments[i]
            seg.update_celerity_diffusivity()
            dt = min(dt, np.min(np.abs(seg.geo['dx'][:] / seg.celerity[:])))
        return self.CFL * dt

    def step(self, theta=0.5):
        """one time step; returns dt"""
        if self._has_junctions:
            return self._step_network(theta)
        dt = self._time_step()
        self.time += dt
        for i in self.calcOrder:                            # forward (discharge) sequence
            seg = self.segments[i]
            seg.read_laterals(self.time, dt)
            seg.Q[0] = seg.read_upstream_Q(self.time)
            seg.solve_discharge(dt, theta)
        self._backward_sequence()                           # backward (depth) sequence
        for i in range(self.num_segments):
            self.segments[i].oldQ[:] = self.segments[i].Q[:]
        return dt

    # ---- network: volume-balance junction condition ----------------------------------------------------------------
    def _prepare_volume_junction(self):
        ds = self.segDownstreamInfo
        self._tribs = [i for i in self.calcOrder if len(ds[i]) > 0]       # channels draining into a junction
        self._tidx = {i: k for k, i in enumerate(self._tribs)}
        self._carry = {i: 0.0 for i in self._tribs}       # volume residual carried to the next step (m3)
        self._lat_old = {i: 0.0 for i in self._tribs}
        self._S = {i: self.segments[i].storage() for i in self._tribs}
        for s in self.segments.values():
            s.sensitivities()
        self.junction_residual = 0.0

    def _outlet_dh(self, seg, dQout):
        """change of the outlet depth over the step (affine in Q_J): rating curve or prescribed depth"""
        if self.dsbcSpec == 0:
            q, hq = seg.bar_params[-1][1, :], seg.bar_params[-1][0, :]
            Q = seg.Q[-1]
            eps = max(1e-6, 1e-6 * abs(Q))
            return (np.interp(Q + eps, q, hq) - np.interp(Q - eps, q, hq)) / (2 * eps) * dQout
        out = np.zeros_like(dQout)
        out[0] = seg.read_downstream_h(self.time) - seg.h[-1]
        return out

    def _step_network(self, theta=0.5):
        """one time step of a network with the volume-balance junction condition (module docstring); returns dt"""
        ds, us = self.segDownstreamInfo, self.segUpstreamInfo
        if not hasattr(self, '_tribs'):
            self._prepare_volume_junction()
        dt = self._time_step()
        self.time += dt
        tribs, tidx = self._tribs, self._tidx
        m1 = len(tribs) + 1
        # forward: Q = Q_a + Q_in u + Q_J w for every channel, its projections on the storage and on the depth at the
        # upstream end, and the inflow as an affine function of the unknown Q_J (columns 1..m)
        V, P, Qin, lat = {}, {}, {}, {}
        for i in self.calcOrder:
            seg = self.segments[i]
            seg.read_laterals(self.time, dt)
            lat[i] = seg.lateral_volume_rate() if seg.laterals else 0.0
            Vi = seg.cn_operator(dt, theta, len(ds[i]) > 0)
            V[i] = Vi
            a = np.zeros(m1)
            if us[i]:
                for k in us[i]:
                    a[1 + tidx[k]] = 1.0
            else:
                a[0] = seg.read_upstream_Q(self.time)
            Qin[i] = a
            dQa = Vi[:, 0] - seg.Q
            P[i] = (np.matmul(seg.G2, dQa), np.matmul(seg.G2, Vi[:, 1]), np.matmul(seg.G2, Vi[:, 2]), dQa[-1])
        # changes of the depth at the upstream end and of the storage (affine in Q_J), from the outlet up
        h0aff, Saff = {}, {}
        for i in self.calcOrder[::-1]:
            seg = self.segments[i]
            p0, pu, pw, dq_last = P[i]
            base = np.zeros((2, m1))
            base[:, 0] = p0
            base += np.outer(pu, Qin[i])
            if len(ds[i]) > 0:
                base[:, 1 + tidx[i]] += pw
                hd = h0aff[ds[i][0]]
            else:
                dQout = dq_last * np.eye(1, m1, 0)[0] + V[i][-1, 1] * Qin[i]
                hd = self._outlet_dh(seg, dQout)
            h0aff[i] = base[1] + seg.h0H * hd
            Saff[i] = base[0] + seg.gH * hd
        # volume balances of the channels draining into junctions: M Q_J = b
        m = m1 - 1
        Mm = np.empty((m, m))
        b = np.empty(m)
        for i in tribs:
            k = tidx[i]
            seg = self.segments[i]
            row = Saff[i] / dt - theta * Qin[i]
            row[1 + k] += theta
            Mm[k] = row[1:]
            b[k] = -(row[0] - theta * lat[i] - (1 - theta) * (seg.oldQ[0] - seg.oldQ[-1] + self._lat_old[i])
                     + self._carry[i] / dt)
        qa = np.r_[1.0, _solve(Mm, b[:, None])[:, 0]]
        for i in self.calcOrder:
            seg = self.segments[i]
            Vi = V[i]
            seg.Q = Vi[:, 0] + (Qin[i] @ qa) * Vi[:, 1] + ((qa[1 + tidx[i]] * Vi[:, 2]) if len(ds[i]) > 0 else 0.0)
            for j in ds[i]:
                self._jQ[j].set(self.time, seg.Q[-1])
        self._backward_sequence()
        # volume residual of the linearization, carried to the next step
        res = 0.0
        for i in tribs:
            seg = self.segments[i]
            Sn = seg.storage()
            e = (Sn - self._S[i]) - dt * (theta * (seg.Q[0] - seg.Q[-1] + lat[i])
                                          + (1 - theta) * (seg.oldQ[0] - seg.oldQ[-1] + self._lat_old[i])) + self._carry[i]
            self._carry[i] = e
            self._S[i] = Sn
            res = max(res, abs(e) / (dt * max(1.0, abs(seg.Q[0]))))
            self._lat_old[i] = lat[i]
        self.junction_residual = res
        for s in self.segments.values():
            s.sensitivities()
            s.oldQ[:] = s.Q
        return dt

    def _backward_sequence(self):
        """depths from the most downstream channel up: the downstream depth of each channel (rating curve or
        boundary_h at the outlet, junction depth otherwise), then the depth profile; the depth at the upstream end of
        the channel is the downstream depth of the channels draining into it"""
        for i in self.calcOrder[::-1]:
            seg = self.segments[i]
            if i == self.calcOrder[-1] and self.dsbcSpec == 0:
                seg.h[-1] = np.interp(seg.Q[-1], seg.bar_params[-1][1, :], seg.bar_params[-1][0, :])  # normal depth
            else:
                seg.h[-1] = seg.read_downstream_h(self.time)
            seg.solve_depth()
            for j in self.segUpstreamInfo[i]:
                self._jh[j].set(self.time, seg.h[0])

    # ---- output and junction records -------------------------------------------------------------------------------
    def write_state(self, i):
        time_folder = f'{self.caseName}/segment{i}/run/{self.time:.4f}'
        os.makedirs(time_folder, exist_ok=True)
        _write_column(f'{time_folder}/h', self.segments[i].h[:])
        _write_column(f'{time_folder}/Q', self.segments[i].Q[:])

    def reset_junctions(self):
        """junction records: the inflow of every channel that starts at a junction (_jQ) and the downstream depth of
        every channel that ends at a junction (_jh), initialized from the initial state"""
        self._jQ, self._jh = {}, {}
        for i in self.calcOrder:
            for j in self.segDownstreamInfo[i]:
                self._jQ[j] = _JunctionSeries()
            for j in self.segUpstreamInfo[i]:
                self._jh[j] = _JunctionSeries()
                self._jh[j].set(self.time, self.segments[i].h[-1])
        for i in self.calcOrder:
            for j in self.segDownstreamInfo[i]:
                self._jQ[j].add(self.time, self.segments[i].Q[-1])
        for k, seg in self.segments.items():          # the segments read their junction boundaries from these records
            seg._upQ_series = self._jQ.get(k)
            seg._dsh_series = self._jh.get(k)
        self._flush_junction_files()

    def _flush_junction_files(self):
        """mirror the junction records to segment<j>/geo/boundary_Q and boundary_h"""
        for j, s in self._jQ.items():
            np.savetxt(f'{self.caseName}/segment{j}/geo/boundary_Q', s.as_array())
        for j, s in self._jh.items():
            np.savetxt(f'{self.caseName}/segment{j}/geo/boundary_h', s.as_array())

    def trim_junction_records(self, n):
        for s in self._jh.values():
            s.trim(n)
        for s in self._jQ.values():
            s.trim(n)
        self._flush_junction_files()


# ---------------------------------------------------------------------------------------------------------------------
# Single channel
# ---------------------------------------------------------------------------------------------------------------------
class SingleChannel(object):
    """One channel (segment): geometry and tables, RBF collocation operators, discharge and depth solutions."""

    # rows of the uniform-flow tables (bar_params) read in update_celerity_diffusivity: Qbar, Dbar, Cbar, dDbar/dQbar
    _BAR_ROWS = np.array([1, 3, 2, 4])
    _QC_ROWS = np.array([1, 2])     # Qbar, Cbar (depth_profile)
    # rows of the cross-section tables (xsParams): top width; area, top width
    _TOPW_ROW = np.array([1])
    _XS_ROWS = np.array([3, 1])

    def __init__(self, caseName, segmentNo, shpC):
        self.caseName = caseName
        self.segmentNo = segmentNo
        self._upQ_series = None     # set by Network when the upstream boundary is a junction
        self._dsh_series = None     # set by Network when the downstream boundary is a junction
        self._h_boundaries = None   # boundary_h of the outlet, read at the first use
        self.load_geometry()
        self.nodeNo = len(self.geo['nodes'])
        self.xs_param_save()
        self.convey_table_calc()
        self.qbar_cbar_dbar_save()
        self.initialize_conditions()
        self.shpC = shpC
        self.f = np.zeros((self.nodeNo, self.nodeNo + 3))
        self.fx = np.zeros_like(self.f)
        self.fxx = np.zeros_like(self.f)
        self.buildMQ(shpC)
        self._build_fast_tables()

    def _build_fast_tables(self):
        """stacked lookup tables for the vectorized interpolation, and the constant rows of the collocation system"""
        self._bar_tab = np.ascontiguousarray(np.stack(self.bar_params))                                  # (N, 5, 301)
        self._xs_tab = np.ascontiguousarray(np.stack([self.xsParams[k] for k in self.geo['xsInfo']]))    # (N, 4, 301)
        # one table for all lookups of a time step (Qbar, Cbar, Dbar, dDbar/dQbar, top width, area), when the
        # uniform-flow and cross-section tables of every node share their depth grid
        self._tab = None
        if np.array_equal(self._bar_tab[:, 0, :], self._xs_tab[:, 0, :]):
            t = _NodeTable(self._bar_tab[:, 0, :], np.concatenate([self._bar_tab[:, [1, 2, 3, 4], :],
                                                                    self._xs_tab[:, [1, 3], :]], axis=1))
            self._tab = t if t.ok else None
        self._eye = np.eye(self.nodeNo)
        self._dQdx_of = None            # the Q array that dQdx was computed from (solve_depth)
        # rows N..N+2 (moment conditions of the polynomial) and column N of rows 1..N-1 never change
        self.sys[self.nodeNo, :self.nodeNo] = 1
        self.sys[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        self.sys[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2
        self.sys[1:self.nodeNo, self.nodeNo] = 1
        self._x1 = self.geo['nodes'][1:].copy()             # nodes 1..N-1 and their squares
        self._x1sq = self.geo['nodes'][1:] ** 2
        self._bw_st = None              # interval search of the last depth profile (reused for area / top width)
        self._bw_hist = []              # last two depth profiles and their time steps (initial guess of Newton)
        self._dt_cn = None              # time step of the last discharge solution

    def load_geometry(self):
        self.geom_path = f'{self.caseName}/segment{self.segmentNo}/geo/'
        self.geo = {
            'nodes': np.loadtxt(self.geom_path + 'nodes'),
            'slopes': np.loadtxt(self.geom_path + 'slopes'),
            'xsInfo': np.loadtxt(self.geom_path + 'xsInfo', dtype=int),       # cross-section number of every node
            'mannings_n': np.loadtxt(self.geom_path + 'mannings_n')
        }
        gdx = np.zeros_like(self.geo['nodes'])
        gdx[:-1] = (self.geo['nodes'][1:] - self.geo['nodes'][:-1])
        gdx[-1] = gdx[-2]
        self.geo.update({'dx': gdx})

    def initialize_conditions(self):
        run_path = f'{self.caseName}/segment{self.segmentNo}/run/'
        self.Q = np.loadtxt(run_path + '0/Q')
        self.h = np.loadtxt(run_path + '0/h')
        self.Q_boundaries = np.atleast_2d(np.loadtxt(self.geom_path + 'boundary_Q'))
        self.oldQ = self.Q.copy()
        self.laterals = self._load_laterals()
        self.G_new = None               # cumulative lateral flow G at the new and the old time level
        self.G_old = None
        self.dQdx = np.zeros_like(self.Q)
        self.celerity = np.zeros_like(self.Q)
        self.diffusivity = np.zeros_like(self.Q)
        self.Gamma = np.zeros_like(self.Q)
        self.area = np.zeros_like(self.Q)
        self.topW = np.zeros_like(self.Q)
        for i in range(self.nodeNo):
            xs = self.xsParams[self.geo['xsInfo'][i]]
            self.area[i] = np.interp(self.h[i], xs[0, :], xs[3, :])
            self.topW[i] = np.interp(self.h[i], xs[0, :], xs[1, :])
        self.sys = np.zeros((self.nodeNo + 3, self.nodeNo + 3))
        self.rhs = np.zeros(self.nodeNo + 3)

    # ---- RBF collocation matrices ----------------------------------------------------------------------------------
    def buildMQ(self, shpC):
        """f, fx, fxx: MQ basis sqrt(r^2 + c^2), c = shpC * r_av, and the quadratic polynomial (columns N..N+2) and
        their first and second derivatives at the nodes"""
        x = self.geo['nodes']
        N = self.nodeNo
        xdif = x[:, None] - x[None, :]
        r_av = 0
        for d in np.abs(np.diff(x)):
            r_av += d / (N - 1)                             # mean node spacing L / (N - 1)
        r = np.abs(xdif)
        c = shpC * r_av
        self.f[:N, :N] = np.sqrt(r ** 2 + c ** 2)
        self.f[:N, N] = 1
        self.f[:N, N + 1] = x
        self.f[:N, N + 2] = x ** 2
        self.fx[:N, :N] = xdif / self.f[:N, :N]
        self.fx[:N, N + 1] = 1
        self.fx[:N, N + 2] = 2 * x
        self.fxx[:N, :N] = 1 / self.f[:N, :N] - xdif ** 2 / self.f[:N, :N] ** 3
        self.fxx[:N, N + 2] = 2
        SingleChannel._build_operators(self)

    def _build_operators(self):
        """hsys: depth from its derivative at nodes 0..N-2 and its value at the last node (backward sequence);
        fx_invF, fxx_invF: first and second derivative of nodal values"""
        N = self.nodeNo
        S = np.zeros((N + 3, N + 3))
        S[:N - 1, :] = self.fx[:-1, :]
        S[N - 1, :] = self.f[-1, :]
        S[N, :N] = 1
        S[N + 1, :N] = self.geo['nodes'][:]
        S[N + 2, :N] = self.geo['nodes'][:] ** 2
        self.hsys = np.matmul(self.f, np.linalg.pinv(S))
        self.invF = np.linalg.pinv(self.f)
        self.fx_invF = np.matmul(self.fx, self.invF)
        self.fxx_invF = np.matmul(self.fxx, self.invF)

    # ---- lateral flows ---------------------------------------------------------------------------------------------
    def _load_laterals(self):
        """point lateral flows of the channel from geo/laterals (csv: names, x, q_file); q_file in geo/lateralDatas
        holds [time (s), Q_l (m3/s)], negative for an outflow"""
        laterals = []
        cfg_path = os.path.join(self.geom_path, 'laterals')
        if os.path.exists(cfg_path):
            with open(cfg_path, 'r') as f:
                for row in csv.DictReader(f):
                    qfile = row['q_file'].strip()
                    if not os.path.isabs(qfile):
                        qfile = os.path.join(self.geom_path, 'lateralDatas', qfile)
                    laterals.append({'name': row.get('names', 'lat'), 'x': float(row['x']),
                                     'ts': np.atleast_2d(np.loadtxt(qfile))})
        return laterals

    def read_laterals(self, time, dt):
        """cumulative lateral flow G at the old and the new time level of the step"""
        if self.laterals:
            self.G_old = self.lateral_G(time - dt) if self.G_new is None else self.G_new
            self.G_new = self.lateral_G(time)

    def lateral_G(self, time):
        """G(x_i, t) = sum_k Q_l,k(t) W_k(x_i) at the nodes (m3/s)"""
        G = np.zeros(self.nodeNo)
        for L in self.laterals:
            Ql = np.interp(time, L['ts'][:, 0], L['ts'][:, 1])
            if abs(Ql) < 1e-12:
                continue
            if 'W' not in L:
                L['W'] = self._lateral_step(L['x'])
            G += Ql * L['W']
        return G

    def _lateral_step(self, x0):
        """W(x_i) of a point lateral flow at x0: 0 upstream, 1 downstream, 1/2 at a node on x0, 1 at the last node"""
        x = self.geo['nodes']
        W = (x > x0).astype(float)
        W[np.abs(x - x0) <= 1e-9 * (x[-1] - x[0])] = 0.5
        W[-1] = 1.0
        return W

    def lateral_volume_rate(self):
        """lateral flow into the channel over its whole length (m3/s)"""
        return float(self.G_new[-1] - self.G_new[0]) if self.G_new is not None else 0.0

    # ---- parameters of the discharge equation ----------------------------------------------------------------------
    def update_celerity_diffusivity(self):
        """Gamma, C and D at every node from the current Q and h"""
        if self._dQdx_of is not self.Q:    # otherwise computed from this Q at the end of the previous step
            self.dQdx = np.matmul(self.fx_invF, self.Q)
        if self._tab is not None:
            # interval search of the depth profile of the last backward sequence, while self.h is that profile
            st = self._bw_st if (self._bw_st is not None and self._bw_st[0] is self.h) else self._tab.index(self.h)
            Qbar, difbar, celebar, dd_dq = self._tab.values(st, (0, 2, 1, 3), self.h)
        else:
            Qbar, difbar, celebar, dd_dq = _interp_table(self.h, self._bar_tab, self._BAR_ROWS)
        self.Gamma = self.Q / Qbar                     # Q = Qbar(h) Gamma, h from the backward sequence
        k = Qbar / difbar * dd_dq                      # (Qbar / Dbar) dDbar/dQbar
        self.celerity[:] = celebar / 2 * (self.Gamma * (1 + k) + 1 / self.Gamma * (1 - k))
        self.diffusivity[:] = difbar / np.abs(self.Gamma)

    # ---- discharge: theta scheme -----------------------------------------------------------------------------------
    def _set_system(self, dt, theta):
        """collocation matrix of the step: row 0 inflow; rows 1..N-1 the diffusive wave equation, the last of them the
        outflow row dQ/dt + C dQ/dx = 0 (D = 0); rows N..N+2 the moment conditions (constant)"""
        N = self.nodeNo
        self.sys[0, :] = self.f[0, :]
        c1 = self.celerity[1:]
        d1 = self.diffusivity[1:].copy()
        d1[-1] = 0.0
        x1 = self._x1
        self.sys[1:N, :N] = self.f[1:, :N] + dt * theta * (c1[:, None] * self.fx[1:, :N] - d1[:, None] * self.fxx[1:, :N])
        self.sys[1:N, N + 1] = x1 + dt * theta * c1
        self.sys[1:N, N + 2] = 2 * x1 * dt * theta * c1 - 2 * theta * dt * d1 + self._x1sq
        self._dt_cn = dt

    def solve_discharge(self, dt, theta=0.5):
        """Q at the new time level: the inflow Q[0] (set before the call), the diffusive wave equation at the
        interior nodes and the outflow row at the downstream end"""
        N = self.nodeNo
        if self.laterals:
            # P = Q - G: dP/dt + C dP/dx = D d2P/dx2 - dG/dt
            G0, G1 = self.G_old, self.G_new
            P = self.oldQ - G0
            adv = self.celerity * np.matmul(self.fx_invF, P)
            diff = self.diffusivity * np.matmul(self.fxx_invF, P)
            dG = G1 - G0
            self.rhs[0] = self.Q[0] - G1[0]
            self.rhs[1:N - 1] = P[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) - dG[1:-1]
            self.rhs[N - 1] = P[-1] + (1 - theta) * dt * ((-adv)[-1]) - dG[-1]
        else:
            # dQdx was computed from Q at the end of the previous step, and oldQ is that Q (only Q[0] is new)
            adv = self.celerity * (self.dQdx if self._dQdx_of is self.Q else np.matmul(self.fx_invF, self.oldQ))
            diff = self.diffusivity * np.matmul(self.fxx_invF, self.oldQ)
            self.rhs[0] = self.Q[0]
            self.rhs[1:N - 1] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1])
            self.rhs[N - 1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1])
        self._set_system(dt, theta)
        # minimum-norm solution with the cutoff of np.linalg.pinv (singular values <= 1e-15 of the largest dropped)
        u, sv, vt = np.linalg.svd(self.sys)
        coef = np.matmul(self.rhs, u)
        cut = 1e-15 * sv[0]
        if sv[-1] > cut:
            coef /= sv
        else:
            keep = sv > cut
            coef[keep] /= sv[keep]
            coef[~keep] = 0
        self.Q = np.matmul(self.f, np.matmul(coef, vt))
        if self.laterals:
            self.Q = self.Q + G1

    # ---- depth: zero-inertia momentum equation ---------------------------------------------------------------------
    def depth_profile(self, Q, h_start, h_down, tol=1e-9, maxit=50):
        """Depths from dh/dx = S0 (1 - Gamma |Gamma|), Gamma = Q / Qbar(h), with the depth h_down at the downstream
        end, on the collocation system hsys:
               h = hsys [R(h); h_down],   R_i = S0_i (1 - Gamma_i |Gamma_i|),  i = 0 .. N-2.
        Newton's method (Jacobian I - [hsys[:, :N-1] diag(dR/dh), 0], dQbar/dh the slope of the Qbar table in the
        current depth interval) with a backtracking line search on the residual, from h_start.
        Returns (h, Gamma, iterations, largest residual in m)."""
        N = self.nodeNo
        S0 = self.geo['slopes'][:-1]
        Hd = self.hsys[:, :N - 1]
        base = self.hsys[:, N - 1] * h_down
        tab = self._tab
        isfinite = math.isfinite

        def residual(h):
            st = tab.index(h) if tab is not None else None
            if st is not None:
                Qbar, Cbar = tab.values(st, (0,), h)[0], None
            else:
                Qbar, Cbar = _interp_table(h, self._bar_tab, self._QC_ROWS)
            G = Q / Qbar
            r = h - (np.matmul(Hd, S0 * (1 - G[:-1] * np.abs(G[:-1]))) + base)
            nr = float(np.abs(r).max())
            return r, (nr if isfinite(nr) else np.inf), G, Qbar, Cbar, st

        h = np.array(h_start, dtype=float)
        r, nr, G, Qbar, Cbar, st = residual(h)
        it = 0
        while nr > tol and it < maxit:
            it += 1
            if st is not None:
                dQbar = tab.slope[tab.node, 0, st[2]]       # dQbar/dh of the interpolated table
            else:
                dQbar = Cbar * _interp_table(h, self._xs_tab, self._TOPW_ROW)[0]
            J = 2 * S0 * np.abs(G[:-1]) * G[:-1] * dQbar[:-1] / Qbar[:-1]      # dR_i/dh_i
            jac = self._eye.copy()
            jac[:, :N - 1] -= Hd * J[None, :]
            step = np.linalg.solve(jac, -r)
            lam = 1.0
            while True:
                hn = h + lam * step
                trial = residual(hn)
                if trial[1] < nr or lam < 1e-8:
                    break
                lam *= 0.5
            h = hn
            r, nr, G, Qbar, Cbar, st = trial
        self._bw_st = st
        return h, G, it, nr

    def solve_depth(self):
        """depth of the backward sequence from the current downstream depth h[-1], then Gamma, wetted area, top
        width and dQ/dx"""
        # initial guess of Newton's method: the last two depth profiles extrapolated linearly in time (saves
        # iterations; the solution does not depend on it), the current depth otherwise
        h_start = self.h
        hist, dt = self._bw_hist, self._dt_cn
        if dt is not None and len(hist) == 2 and hist[1][1] is not None:
            (h0, _), (h1, dt1) = hist
            guess = h1 + (h1 - h0) * (dt / dt1)
            guess[-1] = self.h[-1]
            if np.isfinite(guess).all() and guess.min() > 0:
                h_start = guess
        self.h, G, self.depth_iterations, self.depth_residual = self.depth_profile(self.Q, h_start, self.h[-1])
        hist.append((self.h.copy(), dt))
        if len(hist) > 2:
            del hist[0]
        self._dt_cn = None
        self.Gamma = G
        if self._bw_st is not None:
            self.area[:], self.topW[:] = self._tab.values(self._bw_st, (5, 4), self.h)
        else:
            self.area[:], self.topW[:] = _interp_table(self.h, self._xs_tab, self._XS_ROWS)
        self.dQdx = np.matmul(self.fx_invF, self.Q)
        self._dQdx_of = self.Q

    # ---- volume-balance junction (Network._step_network) -----------------------------------------------------------
    def prepare_linear(self):
        """constant arrays of the volume-balance junction condition"""
        N = self.nodeNo
        x = self.geo['nodes']
        w = np.zeros(N)
        w[:-1] += (x[1:] - x[:-1]) / 2
        w[1:] += (x[1:] - x[:-1]) / 2
        self._trap_w = w                                    # trapezoid weights of the storage integral
        self._S0m = self.geo['slopes'][:-1].copy()
        self._Hd = np.ascontiguousarray(self.hsys[:, :N - 1])
        self._HdT = np.ascontiguousarray(self._Hd.T)
        self._hdcol = self.hsys[:, N - 1].copy()
        self._rhs2 = np.zeros((N, 2))
        self._rhs2[0, 1] = 1.0
        self._R3 = np.zeros((N + 3, 3))
        self._R3[0, 1] = 1.0
        self.G2 = np.zeros((2, N))

    def storage(self):
        """channel storage, integral of the wetted area (trapezoid over the nodes)"""
        return float(self._trap_w @ self.area)

    def cn_operator(self, dt, theta, junction_end):
        """Collocation system of the step (rows as in solve_discharge; at a junction the last row is Q = Q_J),
        decomposed once. Returns V (N x 3): Q_a (zero inflow, zero Q_J), the unit response to the inflow u and the
        unit response to Q_J w (zero if the channel does not end at a junction)."""
        N = self.nodeNo
        R = self._R3
        if self.laterals:
            G0, G1 = self.G_old, self.G_new
            Pold = self.oldQ - G0
            dG = G1 - G0
            adv = self.celerity * np.matmul(self.fx_invF, Pold)
            diff = self.diffusivity * np.matmul(self.fxx_invF, Pold)
            R[0, 0] = -G1[0]
            R[1:N - 1, 0] = Pold[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) - dG[1:-1]
        else:
            adv = self.celerity * (self.dQdx if self._dQdx_of is self.Q else np.matmul(self.fx_invF, self.oldQ))
            diff = self.diffusivity * np.matmul(self.fxx_invF, self.oldQ)
            R[0, 0] = 0.0
            R[1:N - 1, 0] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1])
        self._set_system(dt, theta)
        sysm = self.sys
        if junction_end:
            row = sysm[N - 1].copy()
            sysm[N - 1, :] = self.f[N - 1, :]
            R[N - 1, 0] = -G1[-1] if self.laterals else 0.0
            R[N - 1, 2] = 1.0
        else:
            if self.laterals:
                R[N - 1, 0] = Pold[-1] + (1 - theta) * dt * ((-adv)[-1]) - dG[-1]
            else:
                R[N - 1, 0] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1])
            R[N - 1, 2] = 0.0
        u, sv, vt = np.linalg.svd(sysm)
        if junction_end:
            sysm[N - 1] = row                   # restore the shared matrix
        keep = sv > 1e-15 * sv[0]
        inv = np.where(keep, 1.0 / np.where(keep, sv, 1.0), 0.0)
        V = np.matmul(self.f, np.matmul(vt.T, inv[:, None] * np.matmul(u.T, R)))
        if self.laterals:
            V[:, 0] += G1                       # Q = P + G
        return V

    def sensitivities(self):
        """Linearization of the depth profile at the current state with the Jacobian of the depth system
        (depth_profile): dh = jac^-1 (H_d diag(dR/dQ) dQ + hsys[:, N-1] dh_down). Gives
        dS = G2[0] . dQ + gH dh_down and dh_0 = G2[1] . dQ + h0H dh_down."""
        N = self.nodeNo
        S0, Hd = self._S0m, self._Hd
        tab = self._tab
        st = self._bw_st if self._bw_st is not None else tab.index(self.h)
        jc = st[2]
        Qbar = tab.values(st, (0,), self.h)[0][:-1]
        dQbar = tab.slope[tab.node, 0, jc][:-1]
        dA = tab.slope[tab.node, 5, jc]                     # dA/dh of the area table (top width)
        G = self.Q[:-1] / Qbar
        aG = np.abs(G)
        jac = self._eye.copy()
        jac[:, :N - 1] -= Hd * (2 * S0 * aG * G * dQbar / Qbar)
        rhs = self._rhs2
        rhs[:, 0] = self._trap_w * dA
        z = _solve(jac.T, rhs)
        dRdQ = -2 * S0 * aG / Qbar
        self.G2[:, :-1] = (np.matmul(self._HdT, z) * dRdQ[:, None]).T
        gh = self._hdcol @ z
        self.gH = float(gh[0])
        self.h0H = float(gh[1])

    # ---- cross-section and uniform-flow tables ---------------------------------------------------------------------
    def _read_xs(self, xs_id):
        """cross-section xs<k>: rows (height above the bed, left bank x, right bank x)"""
        return np.loadtxt(f'{self.caseName}/segment{self.segmentNo}/geo/xSecs/xs{xs_id}')

    def xs_param_save(self):
        """xsParams[k]: rows (depth grid, top width, wetted perimeter, area) of cross-section k"""
        unique_xs = sorted(set(self.geo['xsInfo']))
        self.xsParams = [None] * len(unique_xs)
        for xs_id in unique_xs:
            xs = self._read_xs(xs_id)
            max_h = xs[-1][0]
            heyc = np.linspace(0, max_h, N_TABLE + 1)
            hdif = max_h / N_TABLE
            x1_interp = np.interp(0, xs[:, 0], xs[:, 1])
            x2_interp = np.interp(0, xs[:, 0], xs[:, 2])
            area = np.zeros(len(heyc))
            wp = np.zeros(len(heyc))
            tp = np.zeros(len(heyc))
            wp[0] = x2_interp - x1_interp
            tp[0] = wp[0]
            for i in range(1, len(heyc)):
                x1_interp_n = np.interp(heyc[i], xs[:, 0], xs[:, 1])       # banks at the water surface
                x2_interp_n = np.interp(heyc[i], xs[:, 0], xs[:, 2])
                tp[i] = x2_interp_n - x1_interp_n
                wp[i] = wp[i - 1] + np.sqrt((x1_interp_n - x1_interp) ** 2 + hdif ** 2) \
                    + np.sqrt((x2_interp_n - x2_interp) ** 2 + hdif ** 2)
                area[i] = area[i - 1] + tp[i - 1] * hdif \
                    + np.abs(x1_interp_n - x1_interp) * hdif / 2 \
                    + np.abs(x2_interp_n - x2_interp) * hdif / 2
                x1_interp = x1_interp_n
                x2_interp = x2_interp_n
            self.xsParams[xs_id] = np.vstack([heyc, tp, wp, area])

    def convey_table_calc(self):
        """conveyance K(h) of every node. Up to the bank height H: the section of the xs file. Above it, a compound
        section: the trapezoidal main channel (bottom width b, top width tw at H) and overbank beds rising from H at
        the banks to H + dH at the valley walls (lateral slope dH / Lo, Lo = (twcc - tw) / 2 on each side), then
        vertical walls; the overbank roughness is twice the channel roughness (NWM nCC = 2 n)."""
        self.conveyTable = []
        for i in range(self.nodeNo):
            xs_id = self.geo['xsInfo'][i]
            xs = self._read_xs(xs_id)
            max_h = xs[-1][0]
            heyc = np.linspace(0, max_h, N_TABLE + 1)
            H = xs[1, 0]
            dH = xs[2, 0] - xs[1, 0]
            twcc = xs[-1, 2] - xs[-1, 1]
            b = xs[0, 2] - xs[0, 1]
            tw = xs[1, 2] - xs[1, 1]
            Lo = (twcc - tw) / 2
            y = heyc - H
            ramp = y <= dH
            Acc = np.where(ramp, (twcc - tw) * y ** 2 / (2 * dH), (twcc - tw) * (y - dH / 2))
            sl = np.sqrt(1 + (dH / Lo) ** 2)      # length of the sloped overbank bed per unit horizontal width
            wpcc = np.where(ramp, 2 * Lo * y / dH * sl, 2 * Lo * sl + 2 * (y - dH))
            A = (b + tw) * H / 2 + tw * (heyc - H)
            wp = tw
            rcc = Acc / wpcc
            r = A / wp
            k1 = A * r ** (2 / 3) / self.geo['mannings_n'][i]
            k2 = Acc * rcc ** (2 / 3) / self.geo['mannings_n'][i] / 2         # overbank n = 2 n
            cc = np.zeros((len(heyc), 2))
            cc[:, 0] = heyc
            mask = heyc <= H
            A_xs, P_xs = self.xsParams[xs_id][3, mask], self.xsParams[xs_id][2, mask]
            cc[mask, 1] = A_xs * (A_xs / P_xs) ** (2 / 3) / self.geo['mannings_n'][i]
            cc[~mask, 1] = k1[~mask] + k2[~mask]
            self.conveyTable.append(cc)

    def qbar_cbar_dbar_save(self):
        """bar_params[i]: rows (depth grid, Qbar, Cbar, Dbar, dDbar/dQbar) of node i, uniform flow at the depth:
        Qbar = K sqrt(S0), Cbar = dQbar/dh / B, Dbar = Qbar / (2 B S0)"""
        self.bar_params = [None] * self.nodeNo
        for i in range(self.nodeNo):
            xs_id = self.geo['xsInfo'][i]
            max_h = self._read_xs(xs_id)[-1][0]
            heyc = np.linspace(0, max_h, N_TABLE + 1)
            hdif = max_h / N_TABLE
            qbar = np.zeros(len(heyc))
            dqbar_dh = np.zeros(len(heyc))
            cbar = np.zeros(len(heyc))
            dbar = np.zeros(len(heyc))
            ddbar_dqbar = np.zeros(len(heyc))
            qbar[:] = np.sqrt(self.geo['slopes'][i]) * self.conveyTable[i][:, 1]
            # dQbar/dh: fourth-order central differences, second order next to the ends, one-sided at the ends
            dqbar_dh[2:-2] = (-qbar[4:] + 8 * qbar[3:-1] - 8 * qbar[1:-3] + qbar[0:-4]) / 12 / hdif
            dqbar_dh[1] = (qbar[2] - qbar[0]) / 2 / hdif
            dqbar_dh[-2] = (qbar[-1] - qbar[-3]) / 2 / hdif
            dqbar_dh[0] = (qbar[1] - qbar[0]) / hdif
            dqbar_dh[-1] = (qbar[-1] - qbar[-2]) / hdif
            cbar[:] = dqbar_dh[:] / self.xsParams[xs_id][1, :]
            dbar[:] = qbar[:] / 2 / self.xsParams[xs_id][1, :] / self.geo['slopes'][i]
            # dDbar/dQbar along the tabulated relations Dbar(h) and Qbar(h) (B changes with h, so this is not
            # 1 / (2 B S0), and C keeps both terms of Cappelaere 1997, Eq. 14a)
            ddbar_dqbar[:] = np.gradient(dbar, heyc) / np.gradient(qbar, heyc)
            self.bar_params[i] = np.vstack([heyc, qbar, cbar, dbar, ddbar_dqbar])

    # ---- boundary values -------------------------------------------------------------------------------------------
    def read_upstream_Q(self, time):
        if self._upQ_series is not None:          # junction: sum of the outflows of the channels draining into it
            return self._upQ_series.value_at(time)
        return np.interp(time, self.Q_boundaries[:, 0], self.Q_boundaries[:, 1])      # inflow hydrograph

    def read_downstream_h(self, time):
        if self._dsh_series is not None:          # junction: depth at the upstream end of the downstream channel
            return self._dsh_series.value_at(time)
        if self._h_boundaries is None:            # outlet: prescribed depth (boundary_h)
            self._h_boundaries = np.atleast_2d(np.loadtxt(self.geom_path + 'boundary_h'))
        return np.interp(time, self._h_boundaries[:, 0], self._h_boundaries[:, 1])
