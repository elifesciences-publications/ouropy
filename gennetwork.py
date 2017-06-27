# -*- coding: utf-8 -*-
"""
Main class is GenNetwork
Also implements Population and Connection

@author: DanielM
"""

from neuron import h
import random
import numpy as np
import matplotlib.pyplot as plt
import math


class GenNetwork(object):
    """The GenNetwork class organizes populations and connections to a network.
    It uses methods to create gennetwork.Population and Connection classes.
    The GenNetwork provides the logic common in creating a network model and
    keeping track of all the cells, synapses, netcons involved. It also sets up
    measuring vectors for entire populations.

    Attributes
    ----------
    populations - list
        A list of the populations currently present in the network
    connections - list
        A list of connections currently present in the network

    Methods
    -------
    __init__
    mk_Exp2SynConnection
    mk_PerforantPathStimulation
    mk_Sprouting
    current_clamp
    voltage_recording
    set_numpy_seed
    run_network

    Use cases
    ---------
    >>> GenNetwork([GranuleCell, MossyCell], [500,15])
    Create an unconnected network of 500 Granule and 15 Mossy cells.
    """

    def __init__(self, celltypes=None, cellnums=None):
        """Initialize instance empty or with cell populations.
        See gennetwork.Population for detailed implementation of Population.

        Parameters
        ----------
        celltypes - sequence of GenNeuron instances or subclasses thereof
            the threshold value for action potential detection
        cellnums - numeric sequence
            specifies number of neurons per population. same len as celltypes

        Returns
        -------
        self

        Use Cases
        ---------
        >>> GenNetwork([GranuleCell, MossyCell], [500,15])
        Create a network with 500 granule cells and 15 mossy cells.
        """

        self.populations = []
        self.connections = []
        if celltypes is None:
            return

        if cellnums is None:
            return

        for idx, cell_type in enumerate(celltypes):
            self.populations.append(Population(cell_type, cellnums[idx], self))

    def mk_population(self, cell_type, n_cells):
        """Initialize instance empty or with cell populations.
        See gennetwork.Population for detailed implementation of Population.

        Parameters
        ----------
        celltypes - sequence of GenNeuron instances or subclasses thereof
            the threshold value for action potential detection
        cellnums - numeric sequence
            specifies number of neurons per population. same len as celltypes

        Returns
        -------
        self

        """

        if not hasattr(self, 'populations'):
            self.populations = []

        self.populations.append(Population(cell_type, n_cells, self))

    def mk_Exp2SynConnection(self, pre_pop, post_pop, target_pool, target_segs,
                             divergence, tau1, tau2, e, thr, delay, weight):
        """Initialize instance empty or with cell populations.
        See gennetwork.Population for detailed implementation of Population.

        Parameters
        ----------
        celltypes - sequence of GenNeuron instances or subclasses thereof
            the threshold value for action potential detection
        cellnums - numeric sequence
            specifies number of neurons per population. same len as celltypes

        Returns
        -------
        self

        Use Cases
        ---------
        >>> GenNetwork([GranuleCell, MossyCell], [500,15])
        Create a network with 500 granule cells and 15 mossy cells.
        """
        if not hasattr(self, 'connections'):
            self.connections = []

        self.connections.append(Exp2SynConnection(
                                pre_pop, post_pop, target_pool, target_segs,
                                divergence, tau1, tau2, e, thr, delay, weight))

    def mk_PerforantPathStimulation(self, stim, post_pop, n_targets,
                                    target_segs, tau1, tau2, e,
                                    thr, delay, weight):
        if not hasattr(self, 'connections'):
            self.connections = []

        self.connections.append(PerforantPathStimulation(stim,
                                post_pop, n_targets, target_segs,
                                tau1, tau2, e, thr, delay, weight))

    def mk_PerforantPathPoissonStimulation(self,post_pop, t_pattern, spat_pattern, target_segs,
                tau1, tau2, e, weight):
        if not hasattr(self, 'connections'):
            self.connections = []
        
        self.connections.append(PerforantPathPoissonStimulation(post_pop, t_pattern, spat_pattern, target_segs,
                tau1, tau2, e, weight))

    def voltage_recording(self, cell_type):
        rnd_int = random.randint(0, len(self.cells[cell_type]) - 1)
        soma_v_vec, t_vec = self.cells[cell_type][rnd_int]._voltage_recording()
        return soma_v_vec, t_vec

    def set_numpy_seed(self, seed):
        """
        Allows you to set the seed of the numpy.random generator that
        GenNetwork and other classes in gennetwork are using.
        You can also access the random number.
        """
        np.random.seed(seed)
        return np.random.seed

    def run_network(self, tstop=1000, dt=1):
        raise NotImplemented("run_network is not implemented yet")
        h.tstop = tstop
        h.run()


class Population(object):
    """This is the model of a generic population.
    A population is a number of cells of a specific type derived from
    genneuron.GenNeuron. The GenPopulation object keeps track of all
    incoming and outgoing connections. It is recommended to create Populations
    through the GenNetwork.mk_population interface of a network the population
    is part of.

    Attributes
    ----------
    parent_network - gennetwork.GenNetwork or derived instances
        The network the population takes part in
    cell_type - genneuron.GenNeuron class or subclass thereof
        The cell type making up the population
    cells - list of genneuron.GenNeuron instances
        A list of cells that currently exist within the population
    connections - list of Connection objects
        A list of outgoing and incoming connections

    Methods
    -------
    __init__
    make_cells
    get_cell_number
    record_aps
    plot_aps
    write_aps
    current_clamp_rnd
    current_clamp_range
    voltage_recording
    add_connection

    Use cases
    ---------
    >>> nw = GenNetwork()
    >>> nw.mk_population(GranuleCell, 500)
    Create an empty network and create a population of 500 granule cells in the
    network.
    """

    def __init__(self, cell_type=None, n_cells=None, parent_network=None):

        self.parent_network = parent_network
        self.cell_type = cell_type
        self.cells = []
        self.connections = []
        if cell_type and n_cells:
            self.make_cells(cell_type, n_cells)
        self.i = 0

    def make_cells(self, cell_type, n_cells):
        """Create cells of a certain type

        Parameters
        ----------
        cell_type - genneuron.GenNeuron class of subclass thereof
            the type of the cells to be created
        n_cells - numeric
            number of cells to be created

        Returns
        -------
        None

        Use Cases
        ---------
        >>> popul = Population(parent_network = nw)
        >>> popul.make_cells(GranuleCell, 500)
        Create an empty population within nw and then create 500 granule cells
        """

        if hasattr(self, 'cell_type'):
            if self.cell_type != cell_type:
                raise TypeError("cell_type inconsistent with population")
            else:
                self.cell_type = cell_type

        if not hasattr(self, 'cells'):
            self.cells = []

        for x in range(n_cells):
            self.cells.append(cell_type())

        self.cells = np.array(self.cells, dtype=object)

    def get_cell_number(self):
        return len(self.cells)

    def record_aps(self):
        counters = []
        for cell in self.cells:
            counters.append(cell._AP_counter())

        self.ap_counters = counters
        return counters

    def plot_aps(self):
        cells = []
        for x in self.ap_counters:
            cells.append(x[0].as_numpy())

        # Workaround for matplotlib bug. plt.eventplot throws error when first
        # element empty
        if not np.array(cells[0]).any():
            cells[0] = np.array([0], dtype=float)

        plt.eventplot(cells, linewidth=2)

    def write_aps(self, path):
        ap_list = [x[0].as_numpy() for x in self.ap_counters]
        np.savez(path, *ap_list)

    def current_clamp_rnd(self, n_cells, amp=0.3, dur=5, delay=3):
        """

        """
        chosen_cells = np.random.choice(self.cells, n_cells, replace=False)

        for x in chosen_cells:
            for y in delay:
                x._current_clamp_soma(amp=amp, dur=dur, delay=y)

        return chosen_cells

    def current_clamp_range(self, n_cells, amp=0.3, dur=5, delay=3):
        if type(n_cells) == int:
            n_cells = range(n_cells)

        for cell in n_cells:
            self.cells[cell]._current_clamp_soma(amp=amp, dur=dur, delay=delay)

    def voltage_recording(self, cell_type):
        rnd_int = random.randint(0, len(self.cells) - 1)
        soma_v_vec = self.cells[rnd_int]._voltage_recording()
        return soma_v_vec

    def add_connection(self, conn):
        self.connections.append(conn)

    def __iter__(self):
        return self

    def __getitem__(self, item):
        return self.cells[item]
    
    def __str__(self):
        return str(self.get_cell_number) + ' x' + str(self.cell_type) 

    def next(self):
        if self.i < (len(self.cells)):
            i = self.i
            self.i += 1
            return self.cells[i]
        else:
            self.i = 0
            raise StopIteration()

"""
class Connection(object):
    def __str__(self):
        return str(self.pre_pop) + ' to ' + str(self.post_pop)

    def convergence_avg(self):
        pass
"""

class Exp2SynConnection(object):
    """
    This class connects a pre and a post synaptic population with a Exp2Syn
    synapse.
    """

    def __init__(self, pre_pop, post_pop, target_pool, target_segs, divergence,
                 tau1, tau2, e, thr, delay, weight):
        """
        divergence,
                 tau1, tau2, e, g_max, thr, delay, weight, name = "GC->MC"
                 """
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(),dtype=float) / pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) / post_pop.get_cell_number()) * (2*np.pi)
        
        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []

        for idx, curr_cell_pos in enumerate(pre_pop_pos):

            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos,post_cell_pos))
                
            sort_idc = np.argsort(curr_dist)
            closest_cells = sort_idc[0:target_pool]
            picked_cells = np.random.choice(closest_cells, divergence, replace = False)
            pre_cell_target.append(picked_cells)
            for target_cell in picked_cells:
                
                curr_syns = []
                curr_netcons = []

                curr_seg_pool = post_pop[target_cell].get_segs_by_name(target_segs)
                chosen_seg = np.random.choice(curr_seg_pool)
                for seg in chosen_seg:
                    curr_syn = h.Exp2Syn(chosen_seg(0.5))
                    curr_syn.tau1 = tau1
                    curr_syn.tau2 = tau2
                    curr_syn.e = e
                    curr_syns.append(curr_syn)
                    curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v, curr_syn, thr, delay, weight, sec = pre_pop[idx].soma)
                    curr_netcons.append(curr_netcon)
                    netcons.append(curr_netcons)
                    synapses.append(curr_syns)

        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses

#Possibly Depracate
"""class PerforantPathStimulation(object):

    This class connects a pre and a post synaptic population with a Exp2Syn
    synapse.


    def __init__(self, post_pop, n_targets, target_segs,
                 tau1, tau2, e, thr, delay, weight):

        divergence,
                 tau1, tau2, e, g_max, thr, delay, weight, name = "GC->MC"


        synapses = []
        netcons = []
        netstims = []


        if type(n_targets) == int:
            # Select n_targets from post_pop
            target_cells = np.random.choice(post_pop.cells, n_targets, replace = False)
        else:
            target_cells = post_pop.cells[n_targets]

        for curr_cell in target_cells:
            curr_seg_pool = curr_cell.get_segs_by_name(target_segs)
            for seg in curr_seg_pool:
                chosen_seg = np.random.choice(curr_seg_pool)
                curr_syn = h.Exp2Syn(seg(0.5))
                curr_syn.tau1 = tau1
                curr_syn.tau2 = tau2
                curr_syn.e = e
                curr_netstim = h.NetStim()
                curr_netstim.number = 1
                curr_netstim.start = 3
                curr_netcon = h.NetCon(curr_netstim, curr_syn, thr, delay, weight)
                netstims.append(curr_netstim)
                netcons.append(curr_netcon)
                synapses.append(curr_syn)
            
                
        self.netcons = netcons
        self.pre_cell_targets = np.array(target_cells)
        self.synapses = synapses
        self.netstims = netstims
        # Make the synapse"""

class PerforantPathStimulation(object):
    """
    This class connects a pre and a post synaptic population with a Exp2Syn
    synapse.
    """

    def __init__(self, stim, post_pop, n_targets, target_segs,
                 tau1, tau2, e, thr, delay, weight):
        """
        divergence,
                 tau1, tau2, e, g_max, thr, delay, weight, name = "GC->MC"
        """

        self.pre_pop = stim
        self.post_pop = post_pop
        post_pop.add_connection(self)
        synapses = []
        netcons = []
        
        if type(n_targets) == int:
            # Select n_targets from post_pop
            target_cells = np.random.choice(post_pop.cells, n_targets, replace = False)
        else:
            target_cells = post_pop.cells[n_targets]
        
        for curr_cell in target_cells:
            curr_seg_pool = curr_cell.get_segs_by_name(target_segs)
            for seg in curr_seg_pool:
                #chosen_seg = np.random.choice(curr_seg_pool)
                curr_syn = h.Exp2Syn(seg(0.5))
                curr_syn.tau1 = tau1
                curr_syn.tau2 = tau2
                curr_syn.e = e
                curr_netcon = h.NetCon(stim, curr_syn, thr, delay, weight)
                netcons.append(curr_netcon)
                synapses.append(curr_syn)
                
        self.netcons = netcons
        self.pre_cell_targets = np.array(target_cells)
        self.synapses = synapses
        # Make the synapse
        
class PerforantPathPoissonStimulation(object):
    """
    Patterned Perforant Path stimulation as in Yim et al. 2015.
    uses vecevent.mod -> h.VecStim
    """
    def __init__(self, post_pop, t_pattern, spat_pattern, target_segs,
                tau1, tau2, e, weight):

        post_pop.add_connection(self)
        synapses = []
        netcons = []
        
        target_cells = post_pop.cells[spat_pattern]
        
        pattern = np.random.poisson(10,3).cumsum()
        
        self.vecstim = h.VecStim()
        self.pattern_vec = h.Vector(t_pattern)
        self.vecstim.play(self.pattern_vec)
        
        for curr_cell in target_cells:
            curr_seg_pool = curr_cell.get_segs_by_name(target_segs)
            for seg in curr_seg_pool:
                curr_syn = h.Exp2Syn(seg(0.5))
                curr_syn.tau1 = tau1
                curr_syn.tau2 = tau2
                curr_syn.e = e
                curr_netcon = h.NetCon(self.vecstim, curr_syn)
                curr_netcon.weight[0] = weight
                netcons.append(curr_netcon)
                synapses.append(curr_syn)
                for event in pattern:
                    curr_netcon.event(event)
                    
        self.netcons = netcons
        self.pre_cell_targets = np.array(target_cells)
        self.synapses = synapses

"""HELPERS"""
def pos(rad):
    """
    (x,y) position of a point on a circle with axis origin at (0,0)
    and radius 1.
    x = cx + r * cos(rad) -> x = cos(rad)
    y = cy + r * sin(rad) -> y = sin(rad)
    
    Returns a list of tuples that give the point of each radian passed.
    """
    x_arr = list(np.cos(rad))
    y_arr = list(np.sin(rad))
    
    return [(x_arr[idx],y_arr[idx]) for idx in range(len(x_arr))]

def euclidian_dist(p1,p2):
    """ p1 and p2 must both be of len 2 where p1 = (x1,y1); p2 = (x2,y2)"""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    
    
    
    
    