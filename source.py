import numpy as np
import os, csv
import shutil
from collections import deque, defaultdict
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

class Network(object):
    def __init__(self, caseName, g=9.81):
        self.caseName = caseName
        self.clean_run_directories()
        self.g = g
        self.network = np.loadtxt(caseName + '/network', dtype=int)

        with open(caseName+'/input', "r") as file:
            lines = file.readlines()
        self.CFL = float(lines[1][:-1])
        self.RBFtype = int(lines[4][:-1])
        self.shapeCoef = float(lines[7][:-1])
        self.time_integrator = 0
        self.simEndTime = float(lines[10][:-1])
        self.time = 0
        self.printStep = int(lines[13][:-1])
        self.dsbcSpec = int(lines[16][:-1])

        file.close()

        self.num_segments = sum(
            1 for d in os.listdir(caseName) if d.startswith("segment")
            and os.path.isdir(os.path.join(caseName, d)))
        self.segments = self.load_segments()

        if self.num_segments == 1:
            self.calcOrder = [0]
            self.segUpstreamInfo = defaultdict(list)
            self.segDownstreamInfo = defaultdict(list)
            self.downStreams = []
            self.upStreams = []
        else:
            self.queSegments()

            self.uppermostStreams = np.unique(self.network)
            self.downermostStream = np.unique(self.network)
            remove_vals = np.unique(np.atleast_2d(self.network)[:, 1])
            mask = ~np.isin(self.uppermostStreams, remove_vals)
            self.downStreams = self.uppermostStreams[~mask]
            self.uppermostStreams = self.uppermostStreams[mask]
            remove_vals = np.unique(np.atleast_2d(self.network)[:, 0])
            mask = ~np.isin(self.downermostStream, remove_vals)
            self.upStreams = self.downermostStream[~mask]
            self.downermostStream = self.downermostStream[mask]
        self.junction_downs = [k for k, v in self.segUpstreamInfo.items() if v]

        self.reset_junctions()

    def clean_run_directories(self):
        """Remove all time-step directories in 'run', except '0'."""
        for segment_dir in os.listdir(self.caseName):
            segment_path = os.path.join(self.caseName, segment_dir)

            if os.path.isdir(segment_path) and segment_dir.startswith('segment'):
                run_dir = os.path.join(segment_path, 'run')

                if os.path.isdir(run_dir):
                    for time_dir in os.listdir(run_dir):
                        time_dir_path = os.path.join(run_dir, time_dir)
                        if os.path.isdir(time_dir_path) and time_dir != '0':
                            shutil.rmtree(time_dir_path)

        print(f"Completed cleaning up the 'run' directories in {self.caseName}.")

    def load_segments(self):
        """Initialize and store SingleChannel objects for each segment."""
        segments = {}
        for i in range(self.num_segments):
            segments[i] = SingleChannel(self.caseName, i, self.RBFtype, self.shapeCoef)
        return segments

    def queSegments(self):
        self.segUpstreamInfo = defaultdict(list)
        self.segDownstreamInfo = defaultdict(list)
        in_degree = defaultdict(int)  # Store number of upstream segments for each node
        all_nodes = set()

        # Read network connections (assuming self.network stores tuples of (upstream, downstream))
        for upstream, downstream in self.network:
            self.segUpstreamInfo[downstream].append(upstream)
            self.segDownstreamInfo[upstream].append(downstream)
            in_degree[downstream] += 1
            all_nodes.update([upstream, downstream])

        # Find segments with zero in-degree (no upstream dependencies)
        zero_in_degree = deque([node for node in all_nodes if in_degree[node] == 0])

        # Process in topological order
        self.calcOrder = []
        while zero_in_degree:
            segment = zero_in_degree.popleft()
            self.calcOrder.append(segment)

            for downstream in self.segDownstreamInfo[segment]:
                in_degree[downstream] -= 1
                if in_degree[downstream] == 0:
                    zero_in_degree.append(downstream)

        # Convert to numpy array
        self.calcOrder = np.array(self.calcOrder)

    def warmup(self, tol = 1e-4):
        iter = 0
        if self.time_integrator == 0:
            print('Warming up the initial flow depth.')
            res = np.ones(self.num_segments)
            # for i in self.calcOrder:
            #     run_path = self.caseName + '/segment' + str(i) + '/run/'
            #     self.Q = np.loadtxt(run_path + '0/Q_warmup')
            #     ''''''
            #     self.h = np.loadtxt(run_path + '0/h_warmup')

            for j in range(self.num_segments):
                self.segments[j].h = np.zeros_like(self.segments[j].Q)
                for i in range(self.segments[j].nodeNo):
                    self.segments[j].h[i] = np.interp(self.segments[j].Q[i], self.segments[j].bar_params[i][1, :], self.segments[j].bar_params[i][0, :])
                self.segments[j].h_old[:] = self.segments[j].h[:]

            while True:
                dt = 1e8
                for i in range(self.num_segments):
                    # self.segments[i].readLateral(0)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(0)
                    self.segments[i].update_params_cappaleare(self.segDownstreamInfo[i])
                    # self.segments[i].update_params_cappaleare2(self.segDownstreamInfo[i])

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.nanmin(np.abs(dt_arr[:, 0])))

                dt = self.CFL * dt * .1

                iter += 1

                for i in self.calcOrder:
                    # self.segments[i].readLateral(0)
                    self.segments[i].solveSeg_CN_warmup(dt, self.segDownstreamInfo[i])
                    for j in self.segDownstreamInfo[i]:
                        self.set_junctionbc_null(j)

                    """non-iterative, distinct solution for each segment"""
                for i in self.calcOrder:
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)


                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1] and self.dsbcSpec == 0:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1], self.segments[i].bar_params[-1][1, :],self.segments[i].bar_params[-1][0, :])
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h5()
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h_warmup(self.segments[i].h[0], j)


                for i in range(self.num_segments):
                    res[i] = np.max(np.abs((self.segments[i].h_old[:] - self.segments[i].h[:]) / self.segments[i].h[:]) * 100)
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]

                print(f'Residual (%): {np.max(res):.5f}, Iteration: {iter}')
                if (np.max(res) < tol and iter>10) or iter>999:
                    for i in range(self.num_segments):
                        np.savetxt(f"{self.caseName}/segment{i}/run/0/h", self.segments[i].h[:])
                        # np.savetxt(f"{self.caseName}/segment{i}/run/0/Q", self.segments[i].Q[:])
                        dummy = np.atleast_2d(np.loadtxt(f'{self.caseName}/segment{i}/geo/boundary_h'))
                        dummy[0, :] = np.array([[0, self.segments[i].h[-1]]])
                        np.savetxt(f'{self.caseName}/segment{i}/geo/boundary_h', dummy)
                        Qini = np.loadtxt(f"{self.caseName}/segment{i}/run/0/Q")
                        rmse = np.sqrt(np.sum((self.segments[i].Q - Qini) ** 2 / len(Qini)))
                        print(f'Simulation warmed up in {iter} iterations. Deviation in Q of channel {i}: {rmse:.6f}')
                    break

    def warmup2(self, tol = 1e-2):
        iter = 0
        if self.time_integrator == 0:
            print('Warming up the initial flow depth.')
            res = np.ones(self.num_segments)
            # for i in self.calcOrder:
            #     run_path = self.caseName + '/segment' + str(i) + '/run/'
            #     self.Q = np.loadtxt(run_path + '0/Q_warmup')
            #     ''''''
            #     self.h = np.loadtxt(run_path + '0/h_warmup')

            for j in range(self.num_segments):
                self.segments[j].h = np.zeros_like(self.segments[j].Q)
                for i in range(self.segments[j].nodeNo):
                    self.segments[j].h[i] = np.interp(self.segments[j].Q[i], self.segments[j].bar_params[i][1, :],
                                                      self.segments[j].bar_params[i][0, :])
                self.segments[j].h_old[:] = self.segments[j].h[:]

            while True:
                dt = 1e8
                jbc = np.zeros(self.num_segments)
                ddqddx = np.zeros(self.num_segments)
                for i in range(self.num_segments):
                    # self.segments[i].readLateral(0)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(0)
                    self.segments[i].update_params_cappaleare(self.segDownstreamInfo[i])
                    # self.segments[i].update_params_cappaleare2(self.segDownstreamInfo[i])

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.nanmin(np.abs(dt_arr[:, 0])))


                    # ddqddx[i] = np.matmul(self.segments[i].fxx_invF, self.segments[i].Q)[-1]

                dt = self.CFL * dt * .1

                for i in range(self.num_segments):
                    dx = self.segments[i].geo['nodes'][-1] - self.segments[i].geo['nodes'][-2]
                    coef = dx - self.segments[i].cele[-2] * dt
                    jbc[i] = self.segments[i].Q[-2] + coef / (coef + dx - self.segments[i].cele[-2] * dt) * \
                             (self.segments[i].Q[-1] - self.segments[i].Q[-2])
                iter += 1

                for i in self.calcOrder:
                    # self.segments[i].readLateral(0)
                    self.segments[i].solveSeg_CN(dt, self.segDownstreamInfo[i])
                    for j in self.segDownstreamInfo[i]:
                        self.set_junctionbc_null(j)

                    """non-iterative, distinct solution for each segment"""
                for i in self.calcOrder:
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1] and self.dsbcSpec == 0:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],
                                                           self.segments[i].bar_params[-1][1, :],
                                                           self.segments[i].bar_params[-1][0, :])
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h6(dt)
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h_warmup(self.segments[i].h[0], j)

                for i in range(self.num_segments):
                    res[i] = np.max(
                        np.abs((self.segments[i].h_old[:] - self.segments[i].h[:]) / self.segments[i].h[:]) * 100)
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]

                print(f'Residual (%): {np.max(res):.5f}, Iteration: {iter}')
                if (np.max(res) < tol and iter > 10) or iter > 999:
                    for i in range(self.num_segments):
                        np.savetxt(f"{self.caseName}/segment{i}/run/0/h", self.segments[i].h[:])
                        # np.savetxt(f"{self.caseName}/segment{i}/run/0/Q", self.segments[i].Q[:])
                        dummy = np.atleast_2d(np.loadtxt(f'{self.caseName}/segment{i}/geo/boundary_h'))
                        dummy[0, :] = np.array([[0, self.segments[i].h[-1]]])
                        np.savetxt(f'{self.caseName}/segment{i}/geo/boundary_h', dummy)
                        Qini = np.loadtxt(f"{self.caseName}/segment{i}/run/0/Q")
                        rmse = np.sqrt(np.sum((self.segments[i].Q - Qini) ** 2 / len(Qini)))
                        print(f'Simulation warmed up in {iter} iterations. Deviation in Q of channel {i}: {rmse:.6f}')
                    break

    def solve(self):
        iter = 0
        if self.time_integrator == 0:
            omega = .07
            dt_old = 0
            while self.time < self.simEndTime:
                dt = 1e8
                new_xi = np.zeros(self.num_segments)
                for i in range(self.num_segments):
                    # self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)
                    self.segments[i].update_params_cappaleare(self.segDownstreamInfo[i])
                    # self.segments[i].update_params_cappaleare2(self.segDownstreamInfo[i])

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    # dt_arr[:, 1] = self.segments[i].h[:] * self.segments[i].geo['dx'][:] / (self.segments[i].lat[:] + 1e-8) * .2
                    # dt_arr[:, 1] = self.segments[i].area[:] / (self.segments[i].lat[:] + 1e-8) * .2
                    # dt = min(dt, np.min(np.abs(dt_arr[:, 0])), np.min(np.abs(dt_arr[:, 1])))
                    dt = min(dt, np.min(np.abs(dt_arr[:, 0])))
                    # cf = np.max(np.abs(self.segments[i].lat / self.segments[i].latLimiter))
                    # dt /= max(1, cf)

                    # new_xi[i] = self.segments[i].dQdx[-1]

                dt = self.CFL * dt
                iter += 1
                self.time += dt

                jbc = np.zeros(self.num_segments)
                for i in range(self.num_segments):
                    for j in self.segDownstreamInfo[i]:
                        # jbc[i] = self.segments[i].dQdx[-1] + (self.segments[j].cele[0] * dt) / (
                        #         self.segments[i].geo['nodes'][-1] - self.segments[i].geo['nodes'][-2]
                        #         - self.segments[i].cele[-2] * dt
                        #         + self.segments[j].cele[0] * dt) \
                        #          * (self.segments[i].dQdx[-2] - self.segments[j].dQdx[0])
                        for k in self.segUpstreamInfo[j]:
                            if i != k:
                                cele_f = self.segments[i].cele[-1]
                                dqdx_f = self.segments[j].dQdx[0] - self.segments[k].dQdx[-1]
                                # jbc[i] = self.segments[i].dQdx[-1] + (cele_f * dt) / (
                                #         self.segments[i].geo['nodes'][-1] - self.segments[i].geo['nodes'][-2]
                                #         - self.segments[i].cele[-2] * dt
                                #         + cele_f * dt) \
                                #          * (self.segments[i].dQdx[-2] - dqdx_f)
                                # jbc[i] = self.segments[j].cele[0] * self.segments[j].dQdx[0] \
                                #         -self.segments[i].cele[-1] * self.segments[i].dQdx[-1] \
                                #          - self.segments[k].cele[-1] * self.segments[k].dQdx[-1]
                                dx = self.segments[i].geo['nodes'][-1] - self.segments[i].geo['nodes'][-2]
                                coef = dx - self.segments[i].cele[-2] * dt
                                jbc[i] = self.segments[i].Q[-2] + coef / (coef + dx - self.segments[i].cele[-2] * dt) * \
                                                                    (self.segments[i].Q[-1] - self.segments[i].Q[-2])

                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time, dt)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)
                    self.segments[i].solveSeg_CN(dt, self.segDownstreamInfo[i])
                    # for j in self.segDownstreamInfo[i]:
                    #     self.set_junctionbc_null(j)

                    """non-iterative, distinct solution for each segment"""
                # for i in self.calcOrder:
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1] and self.dsbcSpec == 0:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],self.segments[i].bar_params[-1][1, :], self.segments[i].bar_params[-1][0, :])
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h5()
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)

                correction_counter = 0
                res = np.ones(self.num_segments)

                while correction_counter < 1:
                    correction_counter += 1
                    for i in self.calcOrder:
                        for j in self.segDownstreamInfo[i]:
                            self.set_junctionbc_null(j)
                        dtA = (self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt
                        new_xi[i] = omega * -dtA + (1 - omega) * self.segments[i].dQdx[-1]

                        # self.update_celerity(i, -dtA)
                    for i in self.calcOrder:
                        # self.segments[i].readLateral(self.time)
                        self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                        self.segments[i].solveSeg_CN_sec(dt, self.segDownstreamInfo[i], new_xi[i])

                        """non-iterative, distinct solution for each segment"""
                        for j in self.segDownstreamInfo[i]:
                            self.update_junction_Q(self.segments[i].Q[-1], j)

                    for i in self.calcOrder[::-1]:
                        if i == self.calcOrder[-1] and self.dsbcSpec == 0:
                            '''rating curve for the most downstream boundary'''
                            self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1], self.segments[i].bar_params[-1][1, :],self.segments[i].bar_params[-1][0, :])
                        else:
                            self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                        self.segments[i].solveSeg_h5()
                        for j in self.segUpstreamInfo[i]:
                            self.update_junction_h(self.segments[i].h[0], j)

                for i in range(self.num_segments):
                    # self.segments[i].area_old2 = self.segments[i].area_old[-1]
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].area_old[:] = self.segments[i].area[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]

                if iter % self.printStep == 0:
                    print(f"Time: {self.time:.2f}")
                    print(f'------------- Max Q ----- Min Q ------- Max h --- Min h --- Max dQ/dx --------- Min dQ/dx')
                    for i in self.calcOrder:
                        print(f'Channel {i+1}:    {max(self.segments[i].Q):.2f}      {min(self.segments[i].Q):.2f}      {max(self.segments[i].h):.2f}    '
                              f'  {min(self.segments[i].h):.2f}      {max(self.segments[i].dQdx):.6f}      {min(self.segments[i].dQdx):.6f}      {max(self.segments[i].cele):.6f}')
                        # print(self.segments[i].diffu)
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
                        # plt.plot(iter, self.segments[0].dQdx[-1], marker='x', color='k', markersize=3)
                        # plt.plot(iter, self.segments[1].dQdx[-1], marker='x', color='b', markersize=3)
                        # plt.plot(iter, self.segments[2].dQdx[0], marker='x', color='r', markersize=3)
                        # plt.plot(iter, jbc[1], marker='+', color='k', markersize=3)
                        # plt.plot(iter, self.segments[1].Q[-1], marker='+', color='b', markersize=3)
                        # plt.plot(iter, jbc[2], marker='+', color='r', markersize=3)

                    self.clear_boundaryDocs(self.printStep)
                    print(f'----------------------------------------------\n')
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            # plt.show()

    def solve2(self):
        iter = 0
        if self.time_integrator == 0:
            omega = .07
            dt_old = 0

            while self.time < self.simEndTime:
                dt = 1e8

                if iter > 0:
                    pass

                jbc = np.zeros(self.num_segments)
                for i in range(self.num_segments):
                    # self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)
                    self.segments[i].update_params_cappaleare(self.segDownstreamInfo[i])
                    # self.segments[i].update_params_cappaleare2(self.segDownstreamInfo[i])

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    # dt_arr[:, 1] = self.segments[i].h[:] * self.segments[i].geo['dx'][:] / (self.segments[i].lat[:] + 1e-8) * .2
                    # dt_arr[:, 1] = self.segments[i].area[:] / (self.segments[i].lat[:] + 1e-8) * .2
                    # dt = min(dt, np.min(np.abs(dt_arr[:, 0])), np.min(np.abs(dt_arr[:, 1])))
                    dt = min(dt, np.min(np.abs(dt_arr[:, 0])))
                    # cf = np.max(np.abs(self.segments[i].lat / self.segments[i].latLimiter))
                    # dt /= max(1, cf)

                for i in range(self.num_segments):
                    dx = self.segments[i].geo['nodes'][-1] - self.segments[i].geo['nodes'][-2]
                    coef = dx - self.segments[i].cele[-2] * dt
                    jbc[i] = self.segments[i].Q[-2] + coef / (coef + dx - self.segments[i].cele[-2] * dt) * \
                             (self.segments[i].Q[-1] - self.segments[i].Q[-2])

                dt = self.CFL * dt
                iter += 1
                self.time += dt

                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time, dt)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)
                    self.segments[i].solveSeg_CN(dt, self.segDownstreamInfo[i])
                    # for j in self.segDownstreamInfo[i]:
                    #     self.set_junctionbc_null(j)

                    """non-iterative, distinct solution for each segment"""
                # for i in self.calcOrder:
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1] and self.dsbcSpec == 0:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],self.segments[i].bar_params[-1][1, :], self.segments[i].bar_params[-1][0, :])
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h6(dt)
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)

                for i in range(self.num_segments):
                    # self.segments[i].area_old2 = self.segments[i].area_old[-1]
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].area_old[:] = self.segments[i].area[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]

                if iter % self.printStep == 0:
                    print(f"Time: {self.time:.2f}")
                    print(f'------------- Max Q ----- Min Q ------- Max h --- Min h --- Max dQ/dx --------- Min dQ/dx')
                    for i in self.calcOrder:
                        print(f'Channel {i+1}:    {max(self.segments[i].Q):.2f}      {min(self.segments[i].Q):.2f}      {max(self.segments[i].h):.2f}    '
                              f'  {min(self.segments[i].h):.2f}      {max(self.segments[i].dQdx):.6f}      {min(self.segments[i].dQdx):.6f}      {max(self.segments[i].cele):.6f}')
                        # print(self.segments[i].diffu)
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
                        # plt.plot(iter, self.segments[0].dQdx[-1], marker='o', color='k')
                        # plt.plot(iter, self.segments[1].dQdx[-1], marker='o', color='b')
                        # plt.plot(iter, self.segments[2].dQdx[0], marker='o', color='r')
                        # plt.plot(iter, self.segments[2].dQdx[0], marker='o', color='r')

                    self.clear_boundaryDocs(self.printStep)
                    print(f'----------------------------------------------\n')
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            # plt.show()

    def solve_ydk(self):
        iter = 0
        if self.time_integrator == 0:
            # plt.ion()
            # dd = np.array([])
            # cc = np.array([])
            # tmtm = np.array([])
            while self.time < self.simEndTime:
                dt = 1e8
                for i in range(self.num_segments):
                    # self.segments[i].update_params(self.diffLim)
                    self.segments[i].update_params_cappaleare(self.segDownstreamInfo[i])
                    # self.segments[i].update_params_cappaleare2()
                    # self.segments[i].update_params_cappaleare_yedek()
                    # self.segments[i].update_params_new()

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.nanmin(np.abs(dt_arr[:, 0])))

                    # A = (self.segments[i].cele[:] * self.segments[i].fx_invF - self.segments[i].diffu[:] * self.segments[i].fxx_invF)
                    # A = (self.segments[i].cele[:] * self.segments[i].fx_invF)
                    # B = self.segments[i].diffu[:] * self.segments[i].fxx_invF
                    # normA = np.sqrt(np.max(np.abs(np.linalg.eigvals(np.matmul(np.transpose(A), A)))))
                    # normB = np.sqrt(np.max(np.abs(np.linalg.eigvals(np.matmul(np.transpose(B), B)))))
                    # dt = min(dt, 1 / normA)

                    # dx_l = (self.segments[i].geo['nodes'][-1] - self.segments[i].geo['nodes'][-2])
                    # normB = np.abs(self.segments[i].cele[-1] / dx_l - self.segments[i].diffu[-1] / dx_l ** 2)
                    # dt = min(dt, 1 / normA, 1 / normB)

                dt = self.CFL * dt
                # dt = 20


                # print(dt * self.segments[i].cele[0]/ self.segments[i].geo['dx'][1])
                # self.time += dt
                iter += 1
                self.time += dt

                # self.segments[0].Q[0] = self.segments[0].read_upstream_Q(self.time)
                # self.segments[1].Q[0] = self.segments[1].read_upstream_Q(self.time)


                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time, dt)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time) + getattr(self.segments[i], 'qlat_us', 0.0)

                    # self.downstream_characteristics_plusminus(i, dt)
                    # self.downstream_characteristics_plus(i, dt)
                    # self.downstream_characteristics_plus2(i, dt)
                    # self.downstream_characteristics_minus(i, dt)
                    # self.downstream_characteristics_minus2(i, dt)

                    # self.segments[i].solveSeg_fwEuler(dt, self.time)
                    # self.segments[i].solveSeg_fwEuler2(dt, self.time)
                    # self.segments[i].solveSeg_fwEuler3(dt, self.time)
                    # self.segments[i].solveSeg_fwEuler4(dt, self.time)
                    # self.segments[i].solveSeg_fwEuler5(dt, self.time)
                    # self.segments[i].solveSeg_fwEuler6(dt, self.time)
                    # self.segments[i].solveSeg_fwEuler7(dt, self.time)
                    # self.segments[i].solveSeg_bwEuler(dt, self.time)
                    # self.segments[i].solveSeg_heun(dt, self.time)
                    self.segments[i].solveSeg_CN_ydk(dt)
                    # self.segments[i].solveSeg_CN_2(dt, self.time) # simultaneous network
                    # self.segments[i].solveSeg_CN3(dt, self.time, self.segDownstreamInfo[i], self.vmoc)
                    # self.segments[i].solveSeg_CN5(dt, self.time, self.segDownstreamInfo[i])
                    # self.segments[i].solveSeg_RKo2(dt, self.time)
                    # self.segments[i].solveSeg_RKo4(dt, self.time)

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)


                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1]:
                        '''normal conditions for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1], self.segments[i].bar_params[-1][1, :],self.segments[i].bar_params[-1][0, :])
                        '''fixed flow depth for the most dwnstream boundary'''
                        # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    # self.segments[i].solveSeg_h(self.time)
                    # self.segments[i].solveSeg_h2(self.time)
                    # self.segments[i].solveSeg_h3(self.time)
                    # self.segments[i].solveSeg_h_normal()
                    self.segments[i].solveSeg_h4()
                    # self.segments[i].solveSeg_h5()
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)


                for i in range(self.num_segments):
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                # self.time += dt

                if iter % self.printStep == 0:
                    print(f"Time: {self.time}")
                    # dd = np.append(dd, self.segments[1].diffu[0])
                    # cc = np.append(cc, self.segments[1].cele[0])
                    # peclet = cc / dd * self.segments[1].geo['nodes'][-1]
                    # tmtm = np.append(tmtm, self.time)

                    # plt.plot(dd/1000)
                    # plt.plot(cc)
                    # plt.plot(tmtm, cc / dd * self.segments[1].geo['nodes'][-1])
                    # plt.show()
                    # plt.pause(.1)
                    # plt.cla()
                    for i in self.calcOrder:
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            # np.savetxt(f'{self.caseName}/segment{1}/peclet', np.transpose(np.array((tmtm, peclet))))

    def solve_capp2step(self):
        iter = 0
        if self.time_integrator == 0:
            while self.time < self.simEndTime:
                dt = 1e8
                for i in range(self.num_segments):
                    self.segments[i].update_params_cappaleare()

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.nanmin(np.abs(dt_arr[:, 0])))

                dt = self.CFL * dt
                iter += 1
                self.time += dt

                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                    self.segments[i].solveSeg_CN_2step_1(dt, self.segments[i].oldQ)
                    self.segments[i].solveSeg_CN_2step_2(dt, self.segments[i].Q)

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)


                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1]:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1], self.segments[i].bar_params[-1][1, :],self.segments[i].bar_params[-1][0, :])
                        '''fixed flow depth for the most dwnstream boundary'''
                        # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h4()
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)


                for i in range(self.num_segments):
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]

                if iter % self.printStep == 0:
                    print(f"Time: {self.time}")
                    for i in self.calcOrder:
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])

    def solve_iterative_secant(self):
        iter = 0
        if self.time_integrator == 0:
            while self.time < self.simEndTime:
                dt = 1e8
                for i in range(self.num_segments):
                    self.segments[i].update_params_cappaleare()

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.min(np.abs(dt_arr[:, 0])))

                dt = self.CFL * dt
                iter += 1
                self.time += dt

                tet = 4
                secant_xj_old = -np.ones(self.num_segments) * .5e-7
                secant_xj = np.zeros(self.num_segments)
                secant_fj_old = np.zeros(self.num_segments)
                secant_fj = np.zeros(self.num_segments)

                # for i in self.calcOrder:
                #     for j in self.segDownstreamInfo[i]:
                #         self.set_junctionbc_null(j)
                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                    self.segments[i].solveSeg_CN7(dt, self.segDownstreamInfo[i], secant_xj_old[i])

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1]:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],
                                                           self.segments[i].bar_params[-1][1, :],
                                                           self.segments[i].bar_params[-1][0, :])
                        '''fixed flow depth for the most dwnstream boundary'''
                        # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h5(self.time)
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)

                    secant_fj_old[i] = tet * secant_xj_old[i] + (1 - tet) * self.segments[i].dQdx[-1] + (self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt

                for i in self.calcOrder:
                    for j in self.segDownstreamInfo[i]:
                        self.set_junctionbc_null(j)
                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                    self.segments[i].solveSeg_CN7_sec(self.segDownstreamInfo[i], secant_xj[i])

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1]:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],
                                                           self.segments[i].bar_params[-1][1, :],
                                                           self.segments[i].bar_params[-1][0, :])
                        '''fixed flow depth for the most dwnstream boundary'''
                        # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h5(self.time)
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)

                    secant_fj[i] = tet * secant_xj[i] + (1 - tet) * self.segments[i].dQdx[-1] + (
                                self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt

                correction_counter = 0
                res = np.ones(self.num_segments)
                dumq = []
                for n in range(self.num_segments):
                    dumq.append(np.zeros(self.segments[n].nodeNo))
                # while correction_counter < 5 and np.max(res) > 1e-8:
                # while correction_counter < 5:
                while np.max(res) > 1e-4:
                    correction_counter += 1
                    # print(np.max(res))
                    for i in self.calcOrder:
                        for j in self.segDownstreamInfo[i]:
                            self.set_junctionbc_null(j)

                        dummy = secant_xj[i]
                        secant_xj[i] = secant_xj[i] - 1 * secant_fj[i] * (secant_xj_old[i] - secant_xj[i]) / (secant_fj_old[i] - secant_fj[i])
                        secant_xj_old[i] = dummy

                    for i in self.calcOrder:
                        self.segments[i].readLateral(self.time)
                        self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                        self.segments[i].solveSeg_CN7_sec(self.segDownstreamInfo[i], secant_xj[i])

                        """non-iterative, distinct solution for each segment"""
                        for j in self.segDownstreamInfo[i]:
                            self.update_junction_Q(self.segments[i].Q[-1], j)

                    for i in self.calcOrder[::-1]:
                        if i == self.calcOrder[-1]:
                            '''rating curve for the most downstream boundary'''
                            self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1], self.segments[i].bar_params[-1][1, :],self.segments[i].bar_params[-1][0, :])
                            '''fixed flow depth for the most dwnstream boundary'''
                            # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                        else:
                            self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                        self.segments[i].solveSeg_h5(self.time)
                        for j in self.segUpstreamInfo[i]:
                            self.update_junction_h(self.segments[i].h[0], j)

                        dummy = secant_fj[i]
                        secant_fj[i] = tet * secant_xj[i] + (1 - tet) * self.segments[i].dQdx[-1] + (
                                self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt
                        secant_fj_old[i] = dummy

                        res[i] = np.max(np.abs((dumq[i] - self.segments[i].Q) / self.segments[i].Q))
                        dumq[i][:] = self.segments[i].Q[:]

                for i in range(self.num_segments):
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].area_old[:] = self.segments[i].area[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]

                if iter % self.printStep == 0:
                    print(f"Time: {self.time}")
                    for i in self.calcOrder:
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])

    def solve_iterative_bisect(self):
        iter = 0
        if self.time_integrator == 0:
            while self.time < self.simEndTime:
                dt = 1e8
                for i in range(self.num_segments):
                    self.segments[i].update_params_cappaleare()

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.min(np.abs(dt_arr[:, 0])))

                dt = self.CFL * dt
                iter += 1
                self.time += dt

                tet = .85
                bisect_a = -np.ones(self.num_segments) * 1e-1
                bisect_b = np.ones(self.num_segments) * 2e-2
                bisect_c = np.zeros(self.num_segments)
                bisect_fa = np.zeros(self.num_segments)
                bisect_fb = np.zeros(self.num_segments)
                bisect_fc = np.zeros(self.num_segments)

                # for i in self.calcOrder:
                #     for j in self.segDownstreamInfo[i]:
                #         self.set_junctionbc_null(j)
                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                    self.segments[i].solveSeg_CN7(dt, self.segDownstreamInfo[i], bisect_a[i])

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1]:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],
                                                           self.segments[i].bar_params[-1][1, :],
                                                           self.segments[i].bar_params[-1][0, :])
                        '''fixed flow depth for the most dwnstream boundary'''
                        # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h5(self.time)
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)

                    bisect_fa[i] = tet * bisect_a[i] + (1 - tet) * self.segments[i].dQdx[-1] + (self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt

                for i in self.calcOrder:
                    for j in self.segDownstreamInfo[i]:
                        self.set_junctionbc_null(j)
                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                    self.segments[i].solveSeg_CN7_sec(self.segDownstreamInfo[i], bisect_b[i])

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                for i in self.calcOrder[::-1]:
                    if i == self.calcOrder[-1]:
                        '''rating curve for the most downstream boundary'''
                        self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1],
                                                           self.segments[i].bar_params[-1][1, :],
                                                           self.segments[i].bar_params[-1][0, :])
                        '''fixed flow depth for the most dwnstream boundary'''
                        # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    else:
                        self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                    self.segments[i].solveSeg_h5(self.time)
                    for j in self.segUpstreamInfo[i]:
                        self.update_junction_h(self.segments[i].h[0], j)

                    bisect_fb[i] = tet * bisect_b[i] + (1 - tet) * self.segments[i].dQdx[-1] + (
                                self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt

                correction_counter = 0
                res = np.ones(self.num_segments)
                # while correction_counter < 5 and np.max(res) > 1e-8:
                # while correction_counter < 5:
                while np.max(np.abs(res)) > 1e-8:
                    correction_counter += 1
                    print(np.max(bisect_c))
                    for i in self.calcOrder:
                        for j in self.segDownstreamInfo[i]:
                            self.set_junctionbc_null(j)

                        bisect_c[i] = (bisect_b[i] + bisect_a[i]) / 2

                    for i in self.calcOrder:
                        self.segments[i].readLateral(self.time)
                        self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                        self.segments[i].solveSeg_CN7_sec(self.segDownstreamInfo[i], bisect_c[i])

                        """non-iterative, distinct solution for each segment"""
                        for j in self.segDownstreamInfo[i]:
                            self.update_junction_Q(self.segments[i].Q[-1], j)

                    for i in self.calcOrder[::-1]:
                        if i == self.calcOrder[-1]:
                            '''rating curve for the most downstream boundary'''
                            self.segments[i].h[-1] = np.interp(self.segments[i].Q[-1] / self.segments[i].COR[-1], self.segments[i].bar_params[-1][1, :],self.segments[i].bar_params[-1][0, :])
                            '''fixed flow depth for the most dwnstream boundary'''
                            # self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                        else:
                            self.segments[i].h[-1] = self.segments[i].read_downstream_h(self.time)
                        self.segments[i].solveSeg_h5(self.time)
                        for j in self.segUpstreamInfo[i]:
                            self.update_junction_h(self.segments[i].h[0], j)


                        bisect_fc[i] = tet * bisect_c[i] + (1 - tet) * self.segments[i].dQdx[-1] + (
                                self.segments[i].area[-1] - self.segments[i].area_old[-1]) / dt

                        if bisect_fc[i] * bisect_fa[i] < 0:
                            bisect_b[i] = bisect_c[i]
                        else:
                            bisect_a[i] = bisect_c[i]

                        res[i] = np.abs((bisect_a[i] - bisect_b[i]))/bisect_a[i]

                for i in range(self.num_segments):
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].area_old[:] = self.segments[i].area[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]

                if iter % self.printStep == 0:
                    print(f"Time: {self.time}")
                    for i in self.calcOrder:
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])

    def solve_iterative_juncts_h4(self):
        iter = 0
        if self.time_integrator == 0:
            tet = .7
            omega = .01
            while self.time < self.simEndTime:
                dt = 1e8
                for i in range(self.num_segments):
                    self.segments[i].update_params_cappaleare()

                    dt_arr = np.zeros((self.segments[i].cele.shape[0], 2))
                    dt_arr[:, 0] = self.segments[i].geo['dx'][:] / self.segments[i].cele[:]
                    dt = min(dt, np.min(np.abs(dt_arr[:, 0])))

                dt = self.CFL * dt
                iter += 1
                self.time += dt

                new_xi = np.zeros(self.num_segments)
                dumxi = np.zeros(self.num_segments)

                for i in self.calcOrder:
                    self.segments[i].readLateral(self.time)
                    self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                    self.segments[i].solveSeg_CN7(dt, self.segDownstreamInfo[i], new_xi[i])
                    self.segments[i].solveSeg_h4()

                    """non-iterative, distinct solution for each segment"""
                    for j in self.segDownstreamInfo[i]:
                        self.update_junction_Q(self.segments[i].Q[-1], j)

                    dumxi[i] = new_xi[i]
                    new_xi[i] = - (1 - tet) / tet * self.segments[i].dQdx[-1] + 1 / dt / tet * (
                                self.segments[i].area_old[-1] - self.segments[i].area[-1])
                    new_xi[i] = omega * new_xi[i] + (1 - omega) * dumxi[i]

                correction_counter = 0
                res = np.ones(self.num_segments)

                while correction_counter < 3:
                    correction_counter += 1
                    # print(np.max(new_xi))
                    for i in self.junction_downs:
                        self.set_junctionbc_null(i)
                    for i in self.calcOrder:
                        self.segments[i].readLateral(self.time)
                        # self.segments[i].Q[0] = self.segments[i].read_upstream_Q(self.time)

                        self.segments[i].solveSeg_CN7_sec(self.segDownstreamInfo[i], new_xi[i])
                        self.segments[i].solveSeg_h4()

                        """non-iterative, distinct solution for each segment"""
                        for j in self.segDownstreamInfo[i]:
                            self.update_junction_Q(self.segments[i].Q[-1], j)


                        dumxi[i] = new_xi[i]
                        new_xi[i] = - (1 - tet) / tet * self.segments[i].dQdx[-1] + 1 / dt / tet * (
                                self.segments[i].area_old[-1] - self.segments[i].area[-1])
                        new_xi[i] = omega * new_xi[i] + (1 - omega) * dumxi[i]

                        res[i] = np.max(np.abs((dumxi[i] - new_xi[i]) / new_xi[i]))

                for i in range(self.num_segments):
                    self.segments[i].oldQ[:] = self.segments[i].Q[:]
                    self.segments[i].area_old[:] = self.segments[i].area[:]
                    self.segments[i].h_old[:] = self.segments[i].h[:]
                    # self.segments[i].topW_old[:] = self.segments[i].topW[:]

                if iter % self.printStep == 0:
                    print(f"Time: {self.time}")
                    for i in self.calcOrder:
                        time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
                        os.makedirs(time_folder, exist_ok=True)
                        np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
                        np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])
            time_folder = f"{self.caseName}/segment{i}/run/{self.time:.4f}"
            os.makedirs(time_folder, exist_ok=True)
            np.savetxt(f"{time_folder}/h", self.segments[i].h[:])
            np.savetxt(f"{time_folder}/Q", self.segments[i].Q[:])

    def update_junction_Q(self, Q, segId):
        file_path = f"{self.caseName}/segment{segId}/geo/boundary_Q"
        update_tuple = (self.time, Q)

        try:
            # Try reading the file
            with open(file_path, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            # If the file doesn't exist, initialize lines as empty
            lines = []

        updated_lines = []
        entry_written = False

        for line in lines:
            parts = line.split()
            if float(parts[0]) == update_tuple[0]:  # If time matches, update Q
                parts[1] = str(float(parts[1]) + update_tuple[1])
                entry_written = True
            updated_lines.append(" ".join(parts))

        # If no existing entry was updated, append the new tuple
        if not entry_written:
            updated_lines.append(f"{update_tuple[0]} {update_tuple[1]}")

        # Write back the updated content
        with open(file_path, "w") as file:
            file.write("\n".join(updated_lines) + "\n")

    def set_junctionbc_null(self, segId):
        '''this piece of code sets discharge to 0 if same time appears in the text file. '''
        file_path = f"{self.caseName}/segment{segId}/geo/boundary_Q"
        update_tuple = (self.time, 0)

        try:
            # Read existing lines
            with open(file_path, "r") as file:
                lines = file.readlines()
        except FileNotFoundError:
            # If file doesn't exist, start empty
            lines = []

        updated_lines = []
        found = False

        for line in lines:
            # Split existing line into [time, value]
            parts = line.strip().split()
            if len(parts) >= 2 and float(parts[0]) == update_tuple[0]:
                # Replace line if time matches
                updated_lines.append(f"{update_tuple[0]} {update_tuple[1]}\n")
                found = True
            else:
                updated_lines.append(line)

        if not found:
            # Append new line if no match was found
            updated_lines.append(f"{update_tuple[0]} {update_tuple[1]}\n")

        # Write everything back
        with open(file_path, "w") as file:
            file.writelines(updated_lines)

    def update_junction_h(self, h, segId):
        file_path = f"{self.caseName}/segment{segId}/geo/boundary_h"
        new_line = f"{self.time} {h}\n"

        lines = []
        replaced = False

        try:
            with open(file_path, "r") as file:
                for line in file:
                    if line.strip().startswith(f"{self.time} "):
                        lines.append(new_line)
                        replaced = True
                    else:
                        lines.append(line)
        except FileNotFoundError:
            # File does not exist yet
            pass

        if not replaced:
            lines.append(new_line)

        with open(file_path, "w") as file:
            file.writelines(lines)

    def update_junction_h_warmup(self, h, segId):
        file_path = f"{self.caseName}/segment{segId}/geo/boundary_h"
        new_line = f"{self.time} {h}\n"

        lines = []
        replaced = False

        try:
            with open(file_path, "r") as file:
                for line in file:
                    if line.strip().startswith(f"{self.time} "):
                        lines.append(new_line)
                        replaced = True
                    else:
                        lines.append(line)
        except FileNotFoundError:
            # File does not exist yet
            pass

        if not replaced:
            lines.append(new_line)

        with open(file_path, "w") as file:
            file.writelines(lines)

    def reset_junctions(self):
        for i in self.calcOrder:
            for j in self.segDownstreamInfo[i]:
                fpath = f"{self.caseName}/segment{j}/geo/boundary_Q"
                with open(fpath, "w") as file:
                    pass
                file.close()
            for j in self.segUpstreamInfo[i]:
                fpath = f"{self.caseName}/segment{j}/geo/boundary_h"
                with open(fpath, "w") as file:
                    pass
                file.close()
                self.update_junction_h(self.segments[i].h[-1], j)
        for i in self.calcOrder:
            for j in self.segDownstreamInfo[i]:
                self.update_junction_Q(self.segments[i].Q[-1], j)

    def update_celerity(self, i, dQdx):
        Qbar = np.interp(self.segments[i].h[-1], self.segments[i].bar_params[-1][0, :], self.segments[i].bar_params[-1][1, :])
        celebar = np.interp(self.segments[i].h[-1], self.segments[i].bar_params[-1][0, :], self.segments[i].bar_params[-1][2, :])
        difbar = np.interp(self.segments[i].h[-1], self.segments[i].bar_params[-1][0, :], self.segments[i].bar_params[-1][3, :])
        self.segments[i].COR[-1] = np.sqrt(1 - 2 * difbar / celebar / self.segments[i].Q[-1] * dQdx)
        dd_dq = np.interp(self.segments[i].h[-1], self.segments[i].bar_params[-1][0, :], self.segments[i].bar_params[-1][4, :])
        wp = np.interp(self.segments[i].h[-1], self.segments[i].xsParams[self.segments[i].geo['xsInfo'][-1]][0, :],
                       self.segments[i].xsParams[self.segments[i].geo['xsInfo'][-1]][2, :])
        self.segments[i].hydraulic_depth[-1] = self.segments[i].area[-1] / wp
        self.segments[i].cele[-1] = celebar / 2 * (
                    self.segments[i].COR[-1] * (1 + Qbar / difbar * dd_dq) + 1 / self.segments[i].COR[-1] * (1 - Qbar / difbar * dd_dq))

    def clear_boundaryDocs(self, pS):
        for j in self.upStreams:
            fPath = f'{self.caseName}/segment{j}/geo/boundary_h'
            dummy = np.loadtxt(f'{fPath}')
            dummy2 = dummy[-1, :]
            dummy = dummy[:-pS, :]
            dummy = np.append(dummy, [dummy2], axis=0)
            np.savetxt(f'{fPath}', dummy)
        for j in self.downStreams:
            fPath = f'{self.caseName}/segment{j}/geo/boundary_Q'
            dummy = np.loadtxt(f'{fPath}')
            dummy2 = dummy[-1, :]
            dummy = dummy[:-pS, :]
            dummy = np.append(dummy, [dummy2], axis=0)
            np.savetxt(f'{fPath}', dummy)

class SingleChannel(object):
    """Radial Basis Function Collocation Method for 1D diffusive wave equation."""

    def __init__(self, caseName, segmentNo, rbf_type, shpC, g=9.81):
        self.caseName = caseName
        self.segmentNo = segmentNo
        self.g = g
        self.load_geometry()
        self.nodeNo = len(self.geo['nodes'])
        self.xs_param_save()
        self.convey_table_calc()
        self.qbar_cbar_dbar_save()
        self.initialize_conditions()
        self.RBFtype = rbf_type
        self.shpC = shpC
        self.compute_RBF_matrix()

    def load_geometry(self):
        """Load geometry-related data for the segment."""
        self.geom_path = f"{self.caseName}/segment{self.segmentNo}/geo/"
        self.geo = {
            'nodes': np.loadtxt(self.geom_path + 'nodes'),
            'slopes': np.loadtxt(self.geom_path + 'slopes'),
            'xsInfo': np.loadtxt(self.geom_path + 'xsInfo', dtype=int),
            'mannings_n': np.loadtxt(self.geom_path + 'mannings_n')
        }
        gdx = np.zeros_like(self.geo['nodes'])
        gdx[:-1] = (self.geo['nodes'][1:] - self.geo['nodes'][:-1])
        gdx[-1] = gdx[-2]
        # gdx[0], gdx[-1] = (self.geo['nodes'][1] - self.geo['nodes'][0]) / 2, (self.geo['nodes'][-1] - self.geo['nodes'][-2]) / 2
        self.geo.update({'dx': gdx})

    def initialize_conditions(self):
        run_path = self.caseName + '/segment' + str(self.segmentNo) + '/run/'
        self.Q = np.loadtxt(run_path+'0/Q')
        ''''''
        self.h = np.loadtxt(run_path+'0/h')
        ''''''
        # self.h = np.zeros_like(self.Q)
        # for i in range(self.nodeNo):
        #     self.h[i] = np.interp(self.Q[i], self.bar_params[i][1, :], self.bar_params[i][0, :])
        # # self.h[-1] = np.interp(0, np.atleast_2d(self.boundary_h)[:, 0], np.atleast_2d(self.boundary_h)[:, 1])
        # np.savetxt(run_path+'0/h', self.h)
        ''''''
        self.lat = np.zeros_like(self.Q)
        self.area = np.zeros_like(self.Q)
        self.oldQ = np.zeros_like(self.Q)
        self.dQdx = np.zeros_like(self.Q)
        self.Q_boundaries = np.atleast_2d(np.loadtxt(self.geom_path + 'boundary_Q'))
        self.old_dQdx = np.zeros_like(self.Q)
        self.area_old = np.zeros_like(self.Q)
        self.h_old = np.zeros_like(self.Q)
        self.topW = np.zeros_like(self.Q)
        self.topW_old = np.zeros_like(self.Q)
        self.cele = np.zeros_like(self.Q)
        self.cele_h = np.zeros_like(self.Q)
        self.diffu = np.zeros_like(self.Q)
        self.dhdx = np.zeros_like(self.Q)
        self.COR = np.zeros_like(self.Q)
        self.latLimiter = np.ones_like(self.Q) * 1e-4
        self.Sf = np.zeros_like(self.Q)
        self.sys = np.zeros((self.nodeNo+3, self.nodeNo+3))
        self.rhs = np.zeros(self.nodeNo+3)
        self.I = np.eye(self.nodeNo)
        self.bottom_width = np.zeros_like(self.Q)
        self.hydraulic_depth = np.zeros_like(self.Q)
        for i in range(self.nodeNo):
            xs = np.loadtxt(self.caseName + '/segment' + str(self.segmentNo) + '/geo/xSecs/xs' + str(self.geo['xsInfo'][i]))
            self.bottom_width[i] = xs[0,2] - xs[0,1]
            self.area[i] = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0,:], self.xsParams[self.geo['xsInfo'][i]][3,:])
            self.topW[i] = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0,:], self.xsParams[self.geo['xsInfo'][i]][1,:])
        self.area_old[:] = self.area[:]
        self.topW_old[:] = self.topW[:]
        self.h_old[:] = self.h[:]
        self.oldQ[:] = self.Q[:]
        self.old_dQdx[:] = self.dQdx[:]
        self.area_old2 = self.area_old[-1]

    def compute_RBF_matrix(self):
        # self.f = np.zeros((self.nodeNo, self.nodeNo+3), dtype=np.float64)
        self.f = np.zeros((self.nodeNo, self.nodeNo + 3))
        self.fx = np.zeros_like(self.f)
        self.fxx = np.zeros_like(self.f)

        if self.RBFtype == 0:
            self.buildTPS(self.shpC)
        elif self.RBFtype == 1:
            self.buildMQ(self.shpC)
        elif self.RBFtype == 2:
            self.buildIMQ(self.shpC)

    def buildMQ(self, shpC):
        xdif = np.zeros((self.nodeNo, self.nodeNo))
        r_av = 0
        for i in range(self.nodeNo):
            for j in range(self.nodeNo):
                xdif[i, j] = self.geo['nodes'][i] - self.geo['nodes'][j]
                if i == j + 1:
                    r_av += np.abs(xdif[i, j]) / self.nodeNo
        r = np.abs(xdif)
        c = shpC * r_av

        self.f[:self.nodeNo, :self.nodeNo] = np.sqrt(r[:,:] ** 2 + c ** 2)
        self.f[:self.nodeNo, self.nodeNo] = 1
        self.f[:self.nodeNo, self.nodeNo+1] = self.geo['nodes'][:]
        self.f[:self.nodeNo, self.nodeNo + 2] = self.geo['nodes'][:] ** 2
        print(f'{np.linalg.cond(self.f):.2e}')
        self.fx[:self.nodeNo,:self.nodeNo] = xdif[:,:] / self.f[:self.nodeNo,:self.nodeNo]
        self.fx[:self.nodeNo, self.nodeNo + 1] = 1
        self.fx[:self.nodeNo, self.nodeNo + 2] = 2 * self.geo['nodes'][:]
        self.fxx[:self.nodeNo,:self.nodeNo] = 1 / self.f[:self.nodeNo,:self.nodeNo] - xdif[:,:] ** 2 / self.f[:self.nodeNo,:self.nodeNo] ** 3
        self.fxx[:self.nodeNo, self.nodeNo + 2] = 2

        S = np.zeros((self.nodeNo+3, self.nodeNo+3))
        S[:self.nodeNo-1,:] = self.fx[:-1,:]
        S[self.nodeNo-1,:] = self.f[-1,:]
        S[self.nodeNo, :self.nodeNo] = 1
        S[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        S[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2
        invS = np.linalg.pinv(S)
        self.hsys = np.matmul(self.f, invS)

        self.invF = np.linalg.pinv(self.f)
        self.fx_invF = np.matmul(self.fx, self.invF)
        self.fxx_invF = np.matmul(self.fxx, self.invF)

    def buildIMQ(self, shpC): # 0 for 4 rmin, 1 for .815rav
        xdif = np.zeros((self.nodeNo, self.nodeNo))
        r_av=0
        for i in range(self.nodeNo):
            for j in range(self.nodeNo):
                xdif[i, j] = self.geo['nodes'][i] - self.geo['nodes'][j]
                if i == j + 1:
                    r_av += np.abs(xdif[i, j]) / self.nodeNo
        r = np.abs(xdif)
        c = shpC * r_av
        self.f[:self.nodeNo, :self.nodeNo] = 1 / np.sqrt(r[:,:] ** 2 + c ** 2)
        self.fx[:,:self.nodeNo] = - xdif[:,:] * self.f[:,:self.nodeNo]**3
        self.fxx[:,:self.nodeNo] = 3 * xdif[:,:self.nodeNo]**2 * self.f[:,:self.nodeNo]**5 - self.f[:,:self.nodeNo]**3

        self.f[:self.nodeNo, self.nodeNo] = 1
        self.f[:self.nodeNo, self.nodeNo+1] = self.geo['nodes'][:]
        self.f[:self.nodeNo, self.nodeNo + 2] = self.geo['nodes'][:] ** 2
        print(f'{np.linalg.cond(self.f):.2e}')
        self.fx[:self.nodeNo, self.nodeNo + 1] = 1
        self.fx[:self.nodeNo, self.nodeNo + 2] = 2 * self.geo['nodes'][:]
        self.fxx[:self.nodeNo, self.nodeNo + 2] = 2

        S = np.zeros((self.nodeNo+3, self.nodeNo+3))
        S[:self.nodeNo-1,:] = self.fx[:-1,:]
        S[self.nodeNo-1,:] = self.f[-1,:]
        S[self.nodeNo, :self.nodeNo] = 1
        S[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        S[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2
        invS = np.linalg.pinv(S)
        self.hsys = np.matmul(self.f, invS)

        self.invF = np.linalg.pinv(self.f)
        self.fx_invF = np.matmul(self.fx, self.invF)
        self.fxx_invF = np.matmul(self.fxx, self.invF)

    def buildTPS(self, beta):
        for i in range(self.nodeNo):
            for j in range(self.nodeNo):
                if i != j:
                    r = np.abs(self.geo['nodes'][i] - self.geo['nodes'][j])
                    self.f[i,j] = r ** beta * np.log(r)
                    self.fx[i,j] = (self.geo['nodes'][i] - self.geo['nodes'][j]) * r ** (beta - 2) * (
                                beta * np.log(r) + 1)
                    self.fxx[i,j] = r ** (beta - 2) * (beta * np.log(r) + 1) + (
                            self.geo['nodes'][i] - self.geo['nodes'][j]) ** 2 * r ** (beta - 4) * (
                                              2 * (beta - 1) + beta * (beta - 2) * np.log(r))

        self.f[:self.nodeNo, self.nodeNo] = 1
        self.f[:self.nodeNo, self.nodeNo+1] = self.geo['nodes'][:]
        self.f[:self.nodeNo, self.nodeNo + 2] = self.geo['nodes'][:] ** 2
        print(f'{np.linalg.cond(self.f):.2e}')
        self.fx[:self.nodeNo, self.nodeNo + 1] = 1
        self.fx[:self.nodeNo, self.nodeNo + 2] = 2 * self.geo['nodes'][:]

        self.fxx[:self.nodeNo, self.nodeNo + 2] = 2

        S = np.zeros((self.nodeNo+3, self.nodeNo+3))
        S[:self.nodeNo-1,:] = self.fx[:-1,:]
        S[self.nodeNo-1,:] = self.f[-1,:]
        S[self.nodeNo, :self.nodeNo] = 1
        S[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        S[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2
        invS = np.linalg.pinv(S)
        self.hsys = np.matmul(self.f, invS)

        self.invF = np.linalg.pinv(self.f)
        self.fx_invF = np.matmul(self.fx, self.invF)
        self.fxx_invF = np.matmul(self.fxx, self.invF)

    def readLateral_(self, time):
        """Populate self.lat with physically localized lateral source/sink per unit length (q [m^3/s/m]).

        The method supports two modes:
        - *local*: distribute each lateral flow over a short reach around its station using a compact kernel
        - *upstream_lump*: sum all interior lateral flows and add them to the upstream boundary (Fortran-style)

        Lateral definitions are read once from `{segment}/geo/laterals` if present_baton.
            name,x,sign,spread_m,kernel,q_file
        where:
            x         = station location in the same coordinate as `geo/nodes`
            spread_m  = physical mixing length to distribute the source (default: local dx)
            kernel    = 'tri' (default) or 'gauss'
            q_file    = path to a 2-col file [time, Q (m^3/s)] relative to this segment's `geo/` folder or absolute

        The routine is conservative: sum_i lat[i]*dx[i] equals the sum of the specified lateral discharges at the given time.
        """
        # Ensure laterals are loaded
        if not hasattr(self, 'laterals_loaded'):
            self.laterals = []
            cfg_path = os.path.join(self.geom_path, 'laterals')
            if os.path.exists(cfg_path):
                with open(cfg_path, 'r') as f:
                    rdr = csv.DictReader(f)
                    for row in rdr:
                        try:
                            name = row.get('name', 'lat')
                            x = float(row['x'])
                            spread_m = float(row.get('spread_m', '0')) if row.get('spread_m', '').strip() != '' else 0.0
                            kernel = (row.get('kernel', 'tri') or 'tri').strip().lower()
                            qfile = row['q_file'].strip()
                            if not os.path.isabs(qfile):
                                qfile = os.path.join(self.geom_path, 'lateralDatas', qfile)
                            ts = np.atleast_2d(np.loadtxt(qfile))
                        except Exception as e:
                            raise RuntimeError(f"Failed to parse lateral row {row}: {e}")
                        self.laterals.append(
                            {'name': name, 'x': x, 'spread_m': spread_m, 'kernel': kernel, 'ts': ts, '_prev_Q': 0.0})
            self.laterals_mode = os.environ.get('LATERALS_MODE', 'local')  # 'local' or 'upstream_lump'
            self.laterals_loaded = True

        # Reset source term
        self.lat[:] = 0.0
        qlat_us = 0.0

        if not hasattr(self, 'laterals'):  # nothing to do
            self.qlat_us = 0.0
            return

        x_nodes = self.geo['nodes'].copy()
        dx_nodes = self.geo['dx'].copy()

        for L in getattr(self, 'laterals', []):
            if L is None or 'ts' not in L:
                continue
            # timeseries: first column is time in seconds compatible with boundary files; second is Q in m^3/s
            tarr = L['ts'][:, 0];
            varr = L['ts'][:, 1]
            Ql = np.interp(time, tarr, varr)
            if abs(Ql) < 1e-12:
                continue

            # LOCAL APPLICATION (default)
            x0 = L['x']
            # Find a compact support around x0
            dist = np.abs(x_nodes - x0)
            idx0 = int(np.argmin(dist))
            # the span of the lateral flow is increased so that the term
            # dQ/dx does not increase or decrease too sharply.
            maxCMS_in = np.max(self.Q[:])
            maxCMS_lat = np.max(np.abs(varr))
            spread = float(L.get('spread_m', 0.0))
            # spread = np.abs(maxCMS_lat / maxCMS_in) * self.geo['nodes'][-1]
            # spread = float(L.get('spread_m', 0.0)) or float(dx_nodes[idx0])
            R = max(0.5 * spread, 0.5 * dx_nodes[idx0])
            # R = 0.5 * dx_nodes[idx0]
            # mask = dist <= R
            # if not np.any(mask):
            #     mask[idx0] = True
            mask[0], mask[-1] = False, False
            xs = x_nodes[mask]
            dxs = dx_nodes[mask]
            print(np.sum(dxs))
            if L.get('kernel', 'tri') == 'gauss':
                sigma = max(R / 2.0, 1e-6)
                w = np.exp(-0.5 * ((xs - x0) / sigma) ** 2)
            else:
                # triangular kernel
                w = (R - np.abs(xs - x0)) / R
                w[w < 0.0] = 0.0
            denom = np.sum(w * dxs)
            if denom <= 0.0:
                denom = dx_nodes[idx0]
                w = np.ones(1)
                m2 = np.zeros_like(mask);
                m2[idx0] = True
                mask = m2
                dxs = np.array([denom])
            q_per_len = Ql / denom  # ensures sum(lat*dx) = Ql
            self.lat[mask] += q_per_len * w

    def readLateral(self, time, dt):
        """Populate self.lat with physically localized lateral source/sink per unit length (q [m^3/s/m])."""

        # Ensure laterals are loaded (unchanged)
        if not hasattr(self, 'laterals_loaded'):
            self.laterals = []
            cfg_path = os.path.join(self.geom_path, 'laterals')
            if os.path.exists(cfg_path):
                with open(cfg_path, 'r') as f:
                    rdr = csv.DictReader(f)
                    for row in rdr:
                        try:
                            name = row.get('names', 'lat')
                            x = float(row['x'])
                            qfile = row['q_file'].strip()
                            if not os.path.isabs(qfile):
                                qfile = os.path.join(self.geom_path, 'lateralDatas', qfile)
                            ts = np.atleast_2d(np.loadtxt(qfile))
                        except Exception as e:
                            raise RuntimeError(f"Failed to parse lateral row {row}: {e}")
                        self.laterals.append(
                            {'name': name, 'x': x,
                              'ts': ts, '_prev_Q': 0.0
                            }
                        )
            self.laterals_mode = os.environ.get('LATERALS_MODE', 'local')  # 'local' or 'upstream_lump'
            self.laterals_loaded = True

        # Reset source term
        self.lat[:] = 0.0
        qlat_us = 0.0

        if not hasattr(self, 'laterals'):  # nothing to do
            self.qlat_us = 0.0
            return

        x_nodes = self.geo['nodes'].copy()
        dx_nodes = self.geo['dx'].copy()

        for L in getattr(self, 'laterals', []):
            if L is None or 'ts' not in L:
                continue
            # timeseries: first column is time in seconds compatible with boundary files; second is Q in m^3/s
            tarr = L['ts'][:, 0];
            varr = L['ts'][:, 1]
            Ql = np.interp(time, tarr, varr)
            if abs(Ql) < 1e-12:
                continue

            x0 = L['x']
            dist = np.abs(x_nodes - x0)
            idx0 = np.argmin(dist)

            sigmaC = 2.5e-2

            L_tot = x_nodes[-1] - x_nodes[0]  # total length of active reach
            sigma = sigmaC * L_tot  # base width

            # --- NEW: asymmetric sigmas ---
            skew_factor = 4  # >1 means more spread downstream

            sigma_up = sigma  # upstream width  (x <= x0)
            sigma_dn = sigma * skew_factor  # downstream width (x > x0)

            z = x_nodes - x0
            w = np.empty_like(x_nodes)

            # Upstream side (x <= x0): decays quickly
            w[z <= 0] = np.exp(-0.5 * (z[z <= 0] / sigma_up) ** 2)

            # Downstream side (x > x0): decays slowly → more mass downstream
            w[z > 0] = np.exp(-0.5 * (z[z > 0] / sigma_dn) ** 2)
            w[idx0] = 1
            # if  L['name'] == 'bohemia':
            #     w[idx0+1] = 1

            # Keep peak at 1 at x0 (or nearest grid point)
            w /= w.max()
            w[0], w[-1] = 0, 0

            # Optional: visualize
            # plt.plot(x_nodes, w)
            # plt.axvline(x0, ls='--')
            # plt.show()

            denom = np.sum(w * dx_nodes)
            q_per_len = Ql / denom

            self.lat[:] += q_per_len * w[:]

        # coef = 1e3 # higher
        # self.lat[:] = np.sign(self.lat[:]) * np.minimum(np.abs(self.lat[:]), self.oldQ[:] / self.cele[:] / dt / coef)

    def update_params(self, diffLim):
        # dQdx = np.matmul(self.fx_invF, self.Q)
        # dhdx = np.matmul(self.fx_invF, self.h)
        for i in range(self.nodeNo):
            wetPerim, self.area[i], tp = self.interp_wet_area(i)
            R = self.area[i] / wetPerim
            self.Sf[i] = self.geo['mannings_n'][i] ** 2 / (self.area[i]* R ** (2 / 3)) ** 2 * self.Q[i] ** 2
            '''comment out the line above and in the line below to switch wide rectangular channel assumption'''
            # self.Sf[i] = self.geo['mannings_n'][i] ** 2 / (self.area[i] * self.h[i] ** (2 / 3)) ** 2 * self.Q[i] ** 2
            '''Litrico celerity'''
            # self.cele[i] = max(6, 5 * self.Sf[i] ** .3 * self.Q[i] ** .4 / 3 / self.area[i] ** .4 * self.h[i] ** .4 / self.geo['mannings_n'][i] ** .6) if self.h[i] > 0 else 0
            self.cele[i] = 5 * self.Sf[i] ** .3 * self.Q[i] ** .4 / 3 / self.area[i] ** .4 * self.h[i] ** .4 / self.geo['mannings_n'][i] ** .6 if self.h[i] > 0 else 0
            '''Gasiorowski celerity'''
            # self.cele[i] = 5/3 * 1 / self.geo['mannings_n'][i] * R**(2/3) * self.Sf[i]**(.5)
            self.diffu[i] = min(diffLim, np.abs(self.Q[i]) / tp / 2 / self.Sf[i])
        self.cele[:] = np.mean(self.cele[:])
        self.diffu[:] = np.mean(self.diffu[:])

    def update_params_cappaleare(self, dInfo):
        self.dQdx = np.matmul(self.fx_invF, self.Q)
        self.COR = np.ones_like(self.Q)
        for i in range(self.nodeNo - 1):
            Qbar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][1, :])
            difbar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][3, :])
            celebar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][2, :])
            dd_dq = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][4, :])
            # wp = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0, :], self.xsParams[self.geo['xsInfo'][i]][2, :])
            # self.hydraulic_depth[i] = self.area[i] / wp

            # dumcor = np.sqrt(1 - 2 * difbar / celebar / self.Q[i] * self.dQdx[i])
            # self.COR[i] = dumcor if (np.isnan(dumcor) == False) else self.COR[i]
            arg = 1 - 2 * difbar / max(1e-12, celebar) / max(1e-9, self.Q[i]) * self.dQdx[i]
            self.COR[i] = np.sqrt(np.maximum(1e-8, arg))
            # print(self.COR)
            # self.COR[i] = np.sqrt(1 - 2 * difbar / celebar / self.Q[i] * self.dQdx[i])

            self.cele_h[i] = celebar * self.COR[i]
            self.cele[i] = celebar / 2 * (self.COR[i] * (1 + Qbar / difbar * dd_dq) + 1 / self.COR[i] * (1 - Qbar / difbar * dd_dq))
            self.diffu[i] = difbar / self.COR[i]
        Qbar = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][1, :])
        difbar = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][3, :])
        celebar = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][2, :])
        dd_dq = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][4, :])
        # wp = np.interp(self.h[-1], self.xsParams[self.geo['xsInfo'][-1]][0, :],
        #                self.xsParams[self.geo['xsInfo'][-1]][2, :])
        # self.hydraulic_depth[-1] = self.area[-1] / wp
        for i in dInfo:
            # dumcor = np.sqrt(1 - 2 * difbar / celebar / self.Q[i] * self.dQdx[i])
            # self.COR[i] = dumcor if (np.isnan(dumcor) == False) else self.COR[i]
            arg = 1 - 2 * difbar / max(1e-12, celebar) / max(1e-9, self.Q[-1]) * self.dQdx[-1]
            self.COR[-1] = np.sqrt(np.maximum(1e-8, arg))
            # self.COR[-1] = np.sqrt(1 - 2 * difbar / celebar / self.Q[-1] * self.dQdx[-1])

        self.cele_h[-1] = celebar * self.COR[-1]
        self.cele[-1] = celebar / 2 * (self.COR[-1] * (1 + Qbar / difbar * dd_dq) + 1 / self.COR[-1] * (1 - Qbar / difbar * dd_dq))
        self.diffu[-1] = difbar / self.COR[-1]

    def update_params_cappaleare2(self, dInfo):
        self.dQdx = np.matmul(self.fx_invF, self.Q)
        self.COR = np.ones_like(self.Q)
        for i in range(self.nodeNo - 1):
            Qbar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][1, :])
            difbar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][3, :])
            celebar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][2, :])
            dd_dq = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][4, :])
            wp = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0, :], self.xsParams[self.geo['xsInfo'][i]][2, :])
            self.hydraulic_depth[i] = self.area[i] / wp

            self.COR[i] = 1

            self.cele[i] = celebar / 2 * (self.COR[i] * (1 + Qbar / difbar * dd_dq) + 1 / self.COR[i] * (1 - Qbar / difbar * dd_dq))
            self.diffu[i] = difbar / self.COR[i]
        Qbar = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][1, :])
        difbar = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][3, :])
        celebar = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][2, :])
        dd_dq = np.interp(self.h[-1], self.bar_params[-1][0, :], self.bar_params[-1][4, :])
        wp = np.interp(self.h[-1], self.xsParams[self.geo['xsInfo'][-1]][0, :],
                       self.xsParams[self.geo['xsInfo'][-1]][2, :])
        self.hydraulic_depth[-1] = self.area[-1] / wp
        for i in dInfo:
            # dumcor = np.sqrt(1 - 2 * difbar / celebar / self.Q[i] * self.dQdx[i])
            # self.COR[i] = dumcor if (np.isnan(dumcor) == False) else self.COR[i]
            self.COR[-1] = 1

        self.cele[-1] = celebar / 2 * (self.COR[-1] * (1 + Qbar / difbar * dd_dq) + 1 / self.COR[-1] * (1 - Qbar / difbar * dd_dq))
        self.diffu[-1] = difbar / self.COR[-1]

    def update_params_cappaleare_yedek(self):
        dQdx = np.matmul(self.fx_invF, self.Q)
        for i in range(self.nodeNo):
            wetPerim, self.area[i], tp = self.interp_wet_area(i)
            R = self.area[i] / wetPerim
            self.Sf[i] = self.geo['mannings_n'][i] ** 2 / (self.area[i]* R ** (2 / 3)) ** 2 * self.Q[i] ** 2

            Qbar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][1, :])
            difbar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][3, :])
            celebar = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][2, :])
            dd_dq = np.interp(self.h[i], self.bar_params[i][0, :], self.bar_params[i][4, :])
            COR = np.sqrt(1 - 2 * difbar / celebar / self.Q[i] * dQdx[i])
            # COR = self.Q[i] / Qbar
            # COR = 1 if np.isnan(COR) else COR
            # COR = 1

            self.cele[i] = celebar / 2 * (COR * (1 + Qbar / difbar * dd_dq) + 1 / COR * (1 - Qbar / difbar * dd_dq))
            self.diffu[i] = difbar / COR


        # self.cele[:] = 1
        # self.diffu[:] = .01
        # self.cele[:] = np.mean(self.cele[:])
        # self.diffu[:] = np.mean(self.diffu[:])

    def solveSeg_fwEuler(self, dt, time):

        adv = np.matmul(self.I * (self.cele), self.fx_invF)
        diff = np.matmul(self.I * self.diffu, self.fxx_invF)
        lat = self.cele * self.lat

        # self.Q[-1] = self.Q[-2]
        '''Euler here'''
        self.Q[1:-1] += dt * (np.matmul((-adv + diff), self.Q)[1:-1] + lat[1:-1])
        # self.Q[-1] += dt * (np.matmul((-adv), self.Q)[-1] + lat[-1])


        # self.Q[-1] = self.Q[-2]

        # qg = np.exp(self.geo['nodes'][-1] + self.geo['dx'][- 1] - time)
        # self.Q[-1] = qg - self.geo['dx'][- 1] * np.exp(self.geo['nodes'][-1] - time)

        # self.Q[-1] = self.Q[-2] + self.geo['dx'][- 1] * np.exp(self.geo['nodes'][-1] - time)

    def solveSeg_bwEuler(self, dt, time):
        time += dt

        rhs = np.zeros_like(self.Q)
        rhs[0] = self.read_upstream_Q(time)
        rhs[1:-1] = self.Q[1:-1]
        # rhs[-1] = 0
        rhs[-1] = self.Q[-1]

        sys = np.zeros((self.nodeNo, self.nodeNo))
        sys[0,:] = self.f[0,:]
        sys[1:-1,:] = self.f[1:-1,:] + dt * (np.matmul(np.diag(self.cele), self.fx)[1:-1,:] - np.matmul(np.diag(self.diffu), self.fxx)[1:-1,:])
        # sys[-1, :] = self.fx[-1,:]
        sys[-1, :] = self.f[-1,:] + dt * np.matmul(np.diag(self.cele), self.fx)[-1,:]

        invSys = np.linalg.pinv(sys)
        self.Q = np.matmul(self.f, np.matmul(invSys, rhs))

    def solveSeg_bwEuler2(self, dt, time, dInfo):
        # time += dt

        self.rhs[0] = self.read_upstream_Q(time)
        self.rhs[1:-1] = self.oldQ[1:-1]
        self.rhs[-1] = 0

        self.sys = np.zeros((self.nodeNo, self.nodeNo))
        self.sys[0,:] = self.f[0,:]
        self.sys[1:-1,:] = self.f[1:-1,:] + dt * (np.matmul(np.diag(self.cele), self.fx)[1:-1,:] - np.matmul(np.diag(self.diffu), self.fxx)[1:-1,:])
        self.sys[-1, :] = self.fx[-1,:]

        for i in dInfo:
            self.rhs[-1] = -.05 * (self.area[-1] - self.area_old[-1]) / dt

        invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(invSys, self.rhs))

    def solveSeg_heun(self, dt, time):
        adv = np.matmul(np.diag(self.cele), self.fx_invF)
        diff = np.matmul(np.diag(self.diffu), self.fxx_invF)
        lat = self.cele * self.lat

        '''Euler here'''
        oldQ = np.zeros_like(self.Q)
        oldQ[:] = self.Q[:]
        self.Q[1:] = oldQ[1:] + dt * (np.matmul((-adv + diff), oldQ)[1:] + lat[1:])
        self.Q[-1] = self.Q[-2]

        self.solveSeg_h4()
        self.update_params_cappaleare()

        adv = np.matmul(np.diag(self.cele), self.fx_invF)
        diff = np.matmul(np.diag(self.diffu), self.fxx_invF)

        self.Q[1:] = oldQ[1:] + dt * (np.matmul((-adv + diff), self.Q)[1:] + lat[1:])
        self.Q[-1] = self.Q[-2]

        # rhs = np.zeros_like(self.Q)
        # rhs[0] = self.read_upstream_Q(time + dt)
        # rhs[1:-1] = dum[1:-1]
        # rhs[-1] = 0
        #
        # sys = np.zeros((self.nodeNo, self.nodeNo))
        # sys[0,:] = self.f[0,:]
        # sys[1:-1,:] = self.f[1:-1,:] + dt * (np.matmul(np.diag(self.cele), self.fx)[1:-1,:] - np.matmul(np.diag(self.diffu), self.fxx)[1:-1,:])
        # sys[-1, :] = self.fx[-1,:]
        # invSys = np.linalg.pinv(sys)
        # self.Q = np.matmul(self.f, np.matmul(invSys, rhs))

    def solveSeg_CN_ydk(self, dt, theta = .85):
        adv = np.matmul(np.diag(self.cele), np.matmul(self.fx_invF, self.oldQ))
        diff = np.matmul(np.diag(self.diffu), np.matmul(self.fxx_invF, self.oldQ))
        lat = self.cele * self.lat

        self.rhs[0] = self.Q[0]
        self.rhs[1:-1] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1] + lat[1:-1])
        # self.rhs[-1] = - (1 - theta) * self.dQdx[-1]
        self.rhs[-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1] + lat[-1])
        # coef = (self.oldQ[-1] * 2 / self.area_old[-1] + self.cele[-1]) / 2
        # coef = self.oldQ[-1] * 2 / self.area_old[-1]
        # self.rhs[-1] = self.oldQ[-1] - dt * (1 - theta) * coef * self.dQdx[-1]

        carpim_adv = np.matmul(np.diag(self.cele), self.fx)
        carpim_diffu = np.matmul(np.diag(self.diffu), self.fxx)
        self.sys[0, :] = self.f[0, :]
        self.sys[1:-1, :] = self.f[1:-1, :] + dt * theta * (carpim_adv[1:-1, :] - carpim_diffu[1:-1, :])
        # self.sys[-1, :] = theta * self.fx[-1, :]
        self.sys[-1, :] = self.f[-1, :] + dt * theta * (carpim_adv[-1, :])
        # self.sys[-1, :] = self.f[-1, :] + dt * theta * coef * (self.fx[-1, :])

        self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

    def solveSeg_CN(self, dt, dInfo, theta=.85):
        adv = np.matmul(np.diag(self.cele), np.matmul(self.fx_invF, self.oldQ)) # switch with self.dqdx for speed
        diff = np.matmul(np.diag(self.diffu), np.matmul(self.fxx_invF, self.oldQ))
        self.lat = self.cele * self.lat

        self.rhs[0] = self.Q[0]
        self.rhs[1:self.nodeNo-1] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) + dt * self.lat[1:-1]
        # self.rhs[-1] = self.dQdx[-1]
        self.rhs[self.nodeNo-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1])
        # coef = self.oldQ[-1] * 2 / self.area_old[-1]
        # self.rhs[-1] = self.oldQ[-1] - dt * (1 - theta) * coef * self.dQdx[-1]
        # self.latLimiter = np.abs((1 - theta) * ((-adv + diff)))

        carpim_adv = np.matmul(np.diag(self.cele), self.fx[:,:self.nodeNo])
        carpim_diffu = np.matmul(np.diag(self.diffu), self.fxx[:,:self.nodeNo])
        self.sys[0, :] = self.f[0, :]

        self.sys[1:self.nodeNo-1, :self.nodeNo] = self.f[1:-1, :self.nodeNo] + dt * theta * (carpim_adv[1:-1, :] - carpim_diffu[1:-1, :])
        self.sys[1:self.nodeNo - 1, self.nodeNo] = 1
        self.sys[1:self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][1:self.nodeNo - 1] + dt * theta * self.cele[1:self.nodeNo - 1]
        self.sys[1:self.nodeNo - 1, self.nodeNo + 2] = 2 * self.geo['nodes'][1:self.nodeNo - 1] * dt * theta * self.cele[1:self.nodeNo - 1] \
                                                        - 2 * theta * dt * self.diffu[1:self.nodeNo - 1] \
                                                        + self.geo['nodes'][1:self.nodeNo - 1] ** 2

        # self.sys[self.nodeNo-1, :] = self.fx[self.nodeNo-1, :]

        self.sys[self.nodeNo-1, :self.nodeNo] = self.f[self.nodeNo-1, :self.nodeNo] + dt * theta * (carpim_adv[-1, :])
        self.sys[self.nodeNo - 1, self.nodeNo] = 1
        self.sys[self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][-1] + dt * theta * self.cele[-1]
        self.sys[self.nodeNo - 1, self.nodeNo + 2] =  2 * self.geo['nodes'][-1] * dt * theta * self.cele[-1] \
                                                        + self.geo['nodes'][-1] ** 2

        self.sys[self.nodeNo, :self.nodeNo] = 1
        self.sys[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        self.sys[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2

        for i in dInfo:
            # self.rhs[self.nodeNo-1] = 0
            self.rhs[self.nodeNo-1] = self.dQdx[self.nodeNo-1]
            self.sys[self.nodeNo-1, :] = self.fx[-1, :]


        self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

    def solveSeg_CN_sec(self, dt, dInfo, new_xi, theta = .85):
        self.rhs[0] = self.Q[0]
        for i in dInfo:
            self.rhs[self.nodeNo-1] = new_xi
            # self.sys[self.nodeNo-1, :] = self.fx[-1, :]
        # self.rhs[-1] = -(1 - theta) * self.dQdx[-1] + new_xi
        # self.sys[-1, :] = theta * self.fx[-1, :]

        # self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

        # self.apply_lateral_split(dt)

    def solveSeg_CN_warmup(self, dt, dInfo, theta=.85):
        adv = np.matmul(np.diag(self.cele), np.matmul(self.fx_invF, self.oldQ))
        diff = np.matmul(np.diag(self.diffu), np.matmul(self.fxx_invF, self.oldQ))
        lat = self.cele * self.lat

        self.rhs[0] = self.Q[0]
        self.rhs[1:self.nodeNo-1] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) + lat[1:-1] * dt
        self.rhs[self.nodeNo - 1] = 0
        # self.rhs[self.nodeNo-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1])

        carpim_adv = np.matmul(np.diag(self.cele), self.fx[:,:self.nodeNo])
        carpim_diffu = np.matmul(np.diag(self.diffu), self.fxx[:,:self.nodeNo])
        self.sys[0, :] = self.f[0, :]

        self.sys[1:self.nodeNo-1, :self.nodeNo] = self.f[1:-1, :self.nodeNo] + dt * theta * (carpim_adv[1:-1, :] - carpim_diffu[1:-1, :])
        self.sys[1:self.nodeNo - 1, self.nodeNo] = 1
        self.sys[1:self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][1:self.nodeNo - 1] + dt * theta * self.cele[1:self.nodeNo - 1]
        self.sys[1:self.nodeNo - 1, self.nodeNo + 2] = 2 * self.geo['nodes'][1:self.nodeNo - 1] * dt * theta * self.cele[1:self.nodeNo - 1] \
                                                        - 2 * theta * dt * self.diffu[1:self.nodeNo - 1] \
                                                        + self.geo['nodes'][1:self.nodeNo - 1] ** 2

        self.sys[self.nodeNo - 1, :] = self.fx[-1, :]

        # self.sys[self.nodeNo-1, :self.nodeNo] = self.f[self.nodeNo-1, :self.nodeNo] + dt * theta * (carpim_adv[-1, :])
        # self.sys[self.nodeNo - 1, self.nodeNo] = 1
        # self.sys[self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][-1] + dt * theta * self.cele[-1]
        # self.sys[self.nodeNo - 1, self.nodeNo + 2] =  2 * self.geo['nodes'][-1] * dt * theta * self.cele[-1] \
        #                                                 + self.geo['nodes'][-1] ** 2

        self.sys[self.nodeNo, :self.nodeNo] = 1
        self.sys[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        self.sys[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2

        # for i in dInfo:
        #     self.rhs[self.nodeNo-1] = 0
        #     self.sys[self.nodeNo-1, :] = self.fx[-1, :]

        self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

    def solveSeg_CN2(self, dt, dInfo, jbc, theta=.85):
        adv = np.matmul(np.diag(self.cele), np.matmul(self.fx_invF, self.oldQ)) # switch with self.dqdx for speed
        diff = np.matmul(np.diag(self.diffu), np.matmul(self.fxx_invF, self.oldQ))
        self.lat = self.cele * self.lat

        self.rhs[0] = self.Q[0]
        self.rhs[1:self.nodeNo-1] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) + dt * self.lat[1:-1]
        # self.rhs[-1] = self.dQdx[-1]
        self.rhs[self.nodeNo-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1])
        # coef = self.oldQ[-1] * 2 / self.area_old[-1]
        # self.rhs[-1] = self.oldQ[-1] - dt * (1 - theta) * coef * self.dQdx[-1]
        # self.latLimiter = np.abs((1 - theta) * ((-adv + diff)))

        carpim_adv = np.matmul(np.diag(self.cele), self.fx[:,:self.nodeNo])
        carpim_diffu = np.matmul(np.diag(self.diffu), self.fxx[:,:self.nodeNo])
        self.sys[0, :] = self.f[0, :]

        self.sys[1:self.nodeNo-1, :self.nodeNo] = self.f[1:-1, :self.nodeNo] + dt * theta * (carpim_adv[1:-1, :] - carpim_diffu[1:-1, :])
        self.sys[1:self.nodeNo - 1, self.nodeNo] = 1
        self.sys[1:self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][1:self.nodeNo - 1] + dt * theta * self.cele[1:self.nodeNo - 1]
        self.sys[1:self.nodeNo - 1, self.nodeNo + 2] = 2 * self.geo['nodes'][1:self.nodeNo - 1] * dt * theta * self.cele[1:self.nodeNo - 1] \
                                                        - 2 * theta * dt * self.diffu[1:self.nodeNo - 1] \
                                                        + self.geo['nodes'][1:self.nodeNo - 1] ** 2

        # self.sys[self.nodeNo-1, :] = self.fx[self.nodeNo-1, :]

        self.sys[self.nodeNo-1, :self.nodeNo] = self.f[self.nodeNo-1, :self.nodeNo] + dt * theta * (carpim_adv[-1, :])
        self.sys[self.nodeNo - 1, self.nodeNo] = 1
        self.sys[self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][-1] + dt * theta * self.cele[-1]
        self.sys[self.nodeNo - 1, self.nodeNo + 2] =  2 * self.geo['nodes'][-1] * dt * theta * self.cele[-1] \
                                                        + self.geo['nodes'][-1] ** 2

        self.sys[self.nodeNo, :self.nodeNo] = 1
        self.sys[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        self.sys[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2

        for i in dInfo:
            self.rhs[self.nodeNo-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv + diff)[-1])
            self.sys[self.nodeNo - 1, :self.nodeNo] = self.f[-1, :self.nodeNo] + dt * theta * (
                        carpim_adv[-1, :] - carpim_diffu[-1, :])
            self.sys[self.nodeNo - 1, self.nodeNo] = 1
            self.sys[self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][
                                                           self.nodeNo - 1] + dt * theta * self.cele[
                                                                                             self.nodeNo - 1]
            self.sys[self.nodeNo - 1, self.nodeNo + 2] = 2 * self.geo['nodes'][
                                                               self.nodeNo - 1] * dt * theta * self.cele[
                                                                                                 self.nodeNo - 1] \
                                                           - 2 * theta * dt * self.diffu[self.nodeNo - 1] \
                                                           + self.geo['nodes'][self.nodeNo - 1] ** 2

        self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

    def solveSeg_CN_sec2(self, dt, dInfo, convdiff_flux, theta = .85):
        self.rhs[0] = self.Q[0]
        for i in dInfo:
            self.rhs[self.nodeNo-1] = convdiff_flux
            self.sys[self.nodeNo-1, :] = self.cele[-1] * self.f[-1, :] - self.diffu[-1] * self.fx[-1,:]
            # self.sys[self.nodeNo-1, :] = theta * self.fx[-1, :]
        # self.rhs[-1] = -(1 - theta) * self.dQdx[-1] + new_xi
        # self.sys[-1, :] = theta * self.fx[-1, :]

        # self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

        # self.apply_lateral_split(dt)

    def solveSeg_CN_warmup2(self, dt, dInfo, theta=.85):
        adv = np.matmul(np.diag(self.cele), np.matmul(self.fx_invF, self.oldQ))
        diff = np.matmul(np.diag(self.diffu), np.matmul(self.fxx_invF, self.oldQ))
        lat = self.cele * self.lat

        self.rhs[0] = self.Q[0]
        self.rhs[1:self.nodeNo-1] = self.oldQ[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) + lat[1:-1] * dt
        self.rhs[self.nodeNo-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1])
        # self.rhs[-1] = self.oldQ[-1] + (1 - theta) * dt * ((-adv)[-1] + lat[-1])

        carpim_adv = np.matmul(np.diag(self.cele), self.fx[:,:self.nodeNo])
        carpim_diffu = np.matmul(np.diag(self.diffu), self.fxx[:,:self.nodeNo])
        self.sys[0, :] = self.f[0, :]

        self.sys[1:self.nodeNo-1, :self.nodeNo] = self.f[1:-1, :self.nodeNo] + dt * theta * (carpim_adv[1:-1, :] - carpim_diffu[1:-1, :])
        self.sys[1:self.nodeNo - 1, self.nodeNo] = 1
        self.sys[1:self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][1:self.nodeNo - 1] + dt * theta * self.cele[1:self.nodeNo - 1]
        self.sys[1:self.nodeNo - 1, self.nodeNo + 2] = 2 * self.geo['nodes'][1:self.nodeNo - 1] * dt * theta * self.cele[1:self.nodeNo - 1] \
                                                        - 2 * theta * dt * self.diffu[1:self.nodeNo - 1] \
                                                        + self.geo['nodes'][1:self.nodeNo - 1] ** 2

        self.sys[self.nodeNo-1, :self.nodeNo] = self.f[self.nodeNo-1, :self.nodeNo] + dt * theta * (carpim_adv[-1, :])
        self.sys[self.nodeNo - 1, self.nodeNo] = 1
        self.sys[self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][-1] + dt * theta * self.cele[-1]
        self.sys[self.nodeNo - 1, self.nodeNo + 2] =  2 * self.geo['nodes'][-1] * dt * theta * self.cele[-1] \
                                                        + self.geo['nodes'][-1] ** 2

        self.sys[self.nodeNo, :self.nodeNo] = 1
        self.sys[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        self.sys[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2

        for i in dInfo:
            self.rhs[self.nodeNo-1] = jbc
            self.sys[self.nodeNo-1, :] = self.cele[-1] * self.f[-1, :] - self.diffu[-1] * self.fx[-1,:]

        self.invSys = np.linalg.pinv(self.sys)
        self.Q = np.matmul(self.f, np.matmul(self.invSys, self.rhs))

    def solveSeg_RKo2(self, dt, time):
        dummy = np.zeros_like(self.Q)
        K1 = np.zeros_like(self.Q)
        K2 = np.zeros_like(self.Q)
        dummy[:] = self.Q[:]

        '''K1'''
        adv = np.matmul(self.I * (self.cele), np.matmul(self.fx_invF, self.Q))
        diff = np.matmul(self.I * self.diffu, np.matmul(self.fxx_invF, self.Q))
        lat = self.cele * self.lat
        dumt = time + dt
        K1[1:-1] = dt * ((-adv[1:-1] + diff[1:-1]) + lat[1:-1])
        K1[-1] = dt * ((-adv[-1]) + lat[-1])
        dummy[1:] = self.Q[1:] + K1[1:]
        dummy[0] = self.read_upstream_Q(dumt)


        """K2"""
        adv = np.matmul(self.I * (self.cele), np.matmul(self.fx_invF, dummy))
        diff = np.matmul(self.I * self.diffu, np.matmul(self.fxx_invF, dummy))
        lat = self.cele * self.lat
        K2[1:-1] = dt * ((-adv[1:-1] + diff[1:-1]) + lat[1:-1])
        K2[-1] = dt * ((-adv[-1]) + lat[-1])

        self.Q[1:] += 1 / 2 * (K1[1:] + K2[1:])
        self.Q[0] = dummy[0]

    def solveSeg_RKo4(self, dt, time):
        dummy = np.zeros_like(self.Q)
        K1 = np.zeros_like(self.Q)
        K2 = np.zeros_like(self.Q)
        K3 = np.zeros_like(self.Q)
        K4 = np.zeros_like(self.Q)
        dummy[:] = self.Q[:]

        '''K1'''
        adv = np.matmul(self.I * (self.cele), self.fx_invF)
        diff = np.matmul(self.I * self.diffu, self.fxx_invF)
        lat = self.cele * self.lat
        dumt = time + dt / 2
        K1[1:] = dt * (np.matmul((-adv + diff), self.Q)[1:] + lat[1:])
        dummy[1:] = self.Q[1:] + K1[1:] / 2
        dummy[0] = self.read_upstream_Q(dumt)
        # qg = np.exp(self.geo['nodes'][-1] + self.geo['dx'][- 1] - dumt)


        """K2"""
        adv = np.matmul(self.I * (self.cele), self.fx_invF)
        diff = np.matmul(self.I * self.diffu, self.fxx_invF)
        lat = self.cele * self.lat
        dumt = time + dt / 2
        K2[1:] = dt * (np.matmul((-adv + diff), dummy)[1:] + lat[1:])
        dummy[1:] = self.Q[1:] + K2[1:] / 2
        # dummy[0] = self.read_upstream_Q(dumt)
        # qg = np.exp(self.geo['nodes'][-1] + self.geo['dx'][- 1] - dumt)
        # dummy[-1] = qg - self.geo['dx'][- 1] * np.exp(self.geo['nodes'][-1] - dumt)

        """K3"""
        adv = np.matmul(self.I * (self.cele), self.fx_invF)
        diff = np.matmul(self.I * self.diffu, self.fxx_invF)
        lat = self.cele * self.lat
        dumt = time + dt
        K3[1:] = dt * (np.matmul((-adv + diff), dummy)[1:] + lat[1:])
        dummy[1:] = self.Q[1:] + K3[1:]
        dummy[0] = self.read_upstream_Q(dumt)
        # qg = np.exp(self.geo['nodes'][-1] + self.geo['dx'][- 1] - dumt)

        """K4"""
        adv = np.matmul(self.I * (self.cele), self.fx_invF)
        diff = np.matmul(self.I * self.diffu, self.fxx_invF)
        lat = self.cele * self.lat
        K4[1:] = dt * (np.matmul((-adv + diff), dummy)[1:] + lat[1:])

        self.Q[1:] += 1 / 6 * (K1[1:] + 2 * K2[1:] + 2 * K3[1:] + K4[1:])
        self.Q[0] = dummy[0]

    def unit_test_fx_derivative(self):
        # build a linear Q(x) = a*x + b on the collocation nodes
        # you must have self.x or be able to construct coordinates (use node index if needed)
        try:
            x = np.array(self.x)  # prefer actual x coordinates if stored
        except Exception:
            # fallback: assume uniform nodes on [0, L]
            L = getattr(self, 'length', float(self.nodeNo - 1))
            x = np.linspace(0.0, L, self.nodeNo)

        a = 2.34789
        b = -0.1234
        Q_lin = a * x + b

        # derivative estimated by fx row
        slope_fx = float(np.dot(self.fx[-1, :], Q_lin))

        # FD slope (reference): use central diff on nodes near downstream
        # approximate slope using last two node values (simple)
        slope_fd = (Q_lin[-1] - Q_lin[-2]) / (x[-1] - x[-2])

        print(f"[FX TEST] analytic slope a={a:.6e}")
        print(f"[FX TEST] fx-row estimate  = {slope_fx:.6e}")
        print(f"[FX TEST] simple FD slope = {slope_fd:.6e}")
        print(f"[FX TEST] fx - FD = {slope_fx - slope_fd:.6e}")
        return slope_fx, slope_fd

    def solveSeg_h(self, time):
        '''Calculate new h'''
        RHS = self.dhdx
        RHS[-1] = self.h[-1]

        self.h[:-1] = np.matmul(self.hsys[:-1,:], RHS)
        if np.any(self.h < 0):
            for i in range(self.nodeNo - 1):
                self.h[i] = np.interp(self.Q[i], self.bar_params[i][1, :], self.bar_params[i][0, :])
            print('Depths are adjusted with normal depth. Time: ', time)

    def solveSeg_h4(self):
        '''Calculate new h'''
        for i in range(self.nodeNo):
            self.h[i] = np.interp(self.Q[i] / self.COR[i], self.bar_params[i][1, :], self.bar_params[i][0, :])
            self.area[i] = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0, :],
                                 self.xsParams[self.geo['xsInfo'][i]][3, :])

    def solveSeg_h5(self):
        '''Calculate new h'''
        RHS = np.zeros(self.nodeNo+3)
        RHS[:self.nodeNo-1] = self.geo['slopes'][:-1] * (1 - self.COR[:-1] ** 2)
        RHS[self.nodeNo-1] = self.h[-1]
        self.h = np.matmul(self.hsys, RHS)
        for i in range(self.nodeNo):
            self.area[i] = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0, :],
                                 self.xsParams[self.geo['xsInfo'][i]][3, :])
            self.topW[i] = np.interp(self.h[i], self.xsParams[self.geo['xsInfo'][i]][0, :],
                                 self.xsParams[self.geo['xsInfo'][i]][1, :])
        # self.dhdx = np.matmul(self.fx_invF, self.h)
        # self.d2hdx2 = np.matmul(self.fxx_invF, self.h)

    def solveSeg_h6(self, dt, theta = 0.85):
        '''Calculate new h'''
        adv = np.matmul(np.diag(self.cele_h), np.matmul(self.fx_invF, self.h_old))  # switch with self.dqdx for speed
        diff = np.matmul(np.diag(self.diffu), np.matmul(self.fxx_invF, self.h_old))
        self.lat = self.cele_h * self.lat

        self.rhs[0] = self.geo['slopes'][0] * (1 - self.COR[0] ** 2)
        # self.rhs[0] = self.h_old[0] + (1 - theta) * dt * ((-adv + diff)[0])
        # self.rhs[1:self.nodeNo - 1] = self.h_old[1:-1] + (1 - theta) * dt * ((-adv + diff)[1:-1]) # + dt * self.lat[1:-1]
        # self.rhs[self.nodeNo - 1] = self.h[-1]

        carpim_adv = np.matmul(np.diag(self.cele_h), self.fx[:, :self.nodeNo])
        carpim_diffu = np.matmul(np.diag(self.diffu), self.fxx[:, :self.nodeNo])

        self.sys[0, :] = self.fx[0, :]
        # self.sys[0, :] = self.f[0, :]
        # self.sys[0, :self.nodeNo] = self.f[0, :self.nodeNo] + dt * theta * (
        #             carpim_adv[0, :])
        # self.sys[0, self.nodeNo] = 1
        # self.sys[0, self.nodeNo + 1] = self.geo['nodes'][0] + dt * theta * self.cele_h[0]
        # self.sys[0, self.nodeNo + 2] = 2 * self.geo['nodes'][0] * dt * theta * self.cele_h[0] \
        #                                                 + self.geo['nodes'][0] ** 2

        self.sys[1:self.nodeNo - 1, :self.nodeNo] = self.f[1:-1, :self.nodeNo] + dt * theta * (
                    carpim_adv[1:-1, :] - carpim_diffu[1:-1, :])
        self.sys[1:self.nodeNo - 1, self.nodeNo] = 1
        self.sys[1:self.nodeNo - 1, self.nodeNo + 1] = self.geo['nodes'][1:self.nodeNo - 1] + dt * theta * self.cele_h[
                                                                                                           1:self.nodeNo - 1]
        self.sys[1:self.nodeNo - 1, self.nodeNo + 2] = 2 * self.geo['nodes'][
                                                           1:self.nodeNo - 1] * dt * theta * self.cele_h[
                                                                                             1:self.nodeNo - 1] \
                                                       - 2 * theta * dt * self.diffu[1:self.nodeNo - 1] \
                                                       + self.geo['nodes'][1:self.nodeNo - 1] ** 2

        self.sys[self.nodeNo - 1, :] = self.f[-1, :]

        self.sys[self.nodeNo, :self.nodeNo] = 1
        self.sys[self.nodeNo + 1, :self.nodeNo] = self.geo['nodes'][:]
        self.sys[self.nodeNo + 2, :self.nodeNo] = self.geo['nodes'][:] ** 2

        self.invSysH = np.linalg.pinv(self.sys)
        self.h = np.matmul(self.f, np.matmul(self.invSysH, self.rhs))

    def interp_wet_area(self, i):
        h = self.h[i]
        xsNo = self.geo['xsInfo'][i]

        # Interpolate water surface points
        tp = np.interp(h, self.xsParams[xsNo][0,:], self.xsParams[xsNo][1,:])
        wp = np.interp(h, self.xsParams[xsNo][0,:], self.xsParams[xsNo][2,:])
        area = np.interp(h, self.xsParams[xsNo][0,:], self.xsParams[xsNo][3,:])

        return wp, area, tp

    def xs_param_save(self):
        unique_xs = sorted(set(self.geo['xsInfo']))
        self.xsParams = [None] * len(unique_xs)
        NN = 300

        for xs_id in unique_xs:
            xs = np.loadtxt(self.caseName + '/segment' + str(self.segmentNo) + '/geo/xSecs/xs' + str(xs_id))
            max_h = xs[-1][0]
            heyc = np.linspace(0, max_h, NN+1)
            hdif = max_h / NN

            x1_interp = np.interp(0, xs[:, 0], xs[:, 1])
            x2_interp = np.interp(0, xs[:, 0], xs[:, 2])

            area = np.zeros(len(heyc))
            wp = np.zeros(len(heyc))
            tp = np.zeros(len(heyc))

            wp[0] = x2_interp - x1_interp
            tp[0] = wp[0]
            for i in range(1, len(heyc)):

                # Interpolate water surface points

                x1_interp_n = np.interp(heyc[i], xs[:, 0], xs[:, 1])
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
        self.conveyTable = []
        NN = 300

        for i in range(self.nodeNo):
            xs_id = self.geo['xsInfo'][i]
            xs = np.loadtxt(self.caseName + '/segment' + str(self.segmentNo) + '/geo/xSecs/xs' + str(xs_id))
            max_h = xs[-1][0]
            heyc = np.linspace(0, max_h, NN + 1)

            H = xs[1, 0]
            twcc = xs[-1,2] - xs[-1,1]
            b = xs[0,2] - xs[0,1]
            idx_mid = int(len(xs[:,0])/2)
            tw = xs[idx_mid-1,2] - xs[idx_mid-1,1]
            Acc = (twcc - tw) * (heyc - H)
            A = (b + tw) * H / 2 + tw * (heyc - H)
            wpcc = 2 * (heyc - H) + twcc - tw
            wp = tw
            rcc = Acc / wpcc
            r = A / wp
            k1 = A * r ** (2/3) / self.geo['mannings_n'][i]
            k2 = Acc * rcc ** (2/3) / self.geo['mannings_n'][i] / 2 ## ncc is mostly regular mannings doubled, that's why /2

            cc = np.zeros((len(heyc),2))
            cc[:, 0] = heyc
            mask = heyc <= H
            cc[mask, 1] = self.xsParams[xs_id][3,mask] * (self.xsParams[xs_id][3,mask] / self.xsParams[xs_id][2,mask]) ** (2 / 3) / self.geo['mannings_n'][i]
            cc[~mask, 1] = k1[~mask] + k2[~mask]
            self.conveyTable.append(cc)

    def qbar_cbar_dbar_save(self):
        self.bar_params = [None] * self.nodeNo
        NN = 300

        for i in range(self.nodeNo):
            # if i == 52:
            #     continue
            xs_id = self.geo['xsInfo'][i]
            max_h = np.loadtxt(self.caseName + '/segment' + str(self.segmentNo) + '/geo/xSecs/xs' + str(xs_id))[-1][0]
            heyc = np.linspace(0, max_h, NN + 1)
            hdif = max_h / NN

            qbar = np.zeros(len(heyc))
            dqbar_dh = np.zeros(len(heyc))
            cbar = np.zeros(len(heyc))
            dbar = np.zeros(len(heyc))
            ddbar_dqbar = np.zeros(len(heyc))
            # qbar[:] = self.xsParams[xs_id][3, :] * (self.xsParams[xs_id][3, :] / self.xsParams[xs_id][2, :] ) ** (2/3) * np.sqrt(self.geo['slopes'][i]) / self.geo['mannings_n'][i]
            qbar[:] = np.sqrt(self.geo['slopes'][i]) * self.conveyTable[i][:,1]
            '''dqbardh'''
            # dqbar_dh[1:-1] = (qbar[2:] - qbar[:-2]) / 2 / hdif
            ''''''
            dqbar_dh[2:-2] = (-qbar[4:] + 8 * qbar[3:-1] - 8 * qbar[1:-3] + qbar[0:-4]) / 12 / hdif
            dqbar_dh[1] = (qbar[2] - qbar[0]) / 2 / hdif
            dqbar_dh[-2] = (qbar[-1] - qbar[-3]) / 2 / hdif
            ''''''
            dqbar_dh[0] = (qbar[1] - qbar[0]) / hdif
            dqbar_dh[-1] = (qbar[-1] - qbar[-2]) / hdif
            ''''''
            # dqbar_dh[:] = 5 / 3 * heyc[:] ** (2/3) * self.xsParams[xs_id][1, :] * np.sqrt(self.geo['slopes'][i]) / self.geo['mannings_n'][i]
            cbar[:] = dqbar_dh[:] / self.xsParams[xs_id][1, :]
            dbar[:] = qbar[:] / 2 / self.xsParams[xs_id][1, :] / self.geo['slopes'][i]
            ddbar_dqbar[:] = 1 / 2 / self.xsParams[xs_id][1, :] / self.geo['slopes'][i]

            self.bar_params[i] = np.vstack([heyc, qbar, cbar, dbar, ddbar_dqbar])

        # plt.plot(qbar, cbar)
        # plt.plot(qbar, dbar/10/qbar)
        # plt.show()
        # np.savetxt('../Research/Diffusive Wave RBFCM/Q vs CD/qbar_cdbar', np.vstack([qbar,cbar, dbar]))

    def read_upstream_Q(self, time):
        self.Q_boundaries = np.atleast_2d(np.loadtxt(self.geom_path + 'boundary_Q'))
        # Interpolate Q for the given time
        Q_interp = np.interp(time, self.Q_boundaries[:, 0], self.Q_boundaries[:, 1])

        return Q_interp

    def read_downstream_h(self, time):
        h_boundaries = np.atleast_2d(np.loadtxt(self.geom_path + 'boundary_h'))

        # Interpolate Q for the given time
        h_interp = np.interp(time, h_boundaries[:, 0], h_boundaries[:, 1])

        return h_interp