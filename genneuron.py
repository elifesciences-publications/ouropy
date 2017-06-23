# -*- coding: utf-8 -*-
"""
Implementation of a generic neuron.
A neuron always has a soma and a variable number of dendrites.

Uses
"""

from neuron import h, gui
from gendendrite import GenDendrite

class GenNeuron(object):
    """This is the model of a generic neuron.
    It implements those in a generic way so a class of a specific cell type
    should not need to overwrite the components, but could just call the
    methods to initialize it's parameters.
    """

    def __init__(self):
        """Generic __init__ mostly for debugging"""
        self.soma = None
        self.dendrites = None
    
    def mk_soma(self, name = None, diam = None, L = None):
        """Assignes self.soma a hoc section with dimensions diam and L.
        Uses nrn defaults when None. Name defaults to 'soma'.
        Before mk_soma is called, self.soma = None.
        Use cases:
        self.mk_soma()
        - self.soma becomes section with default values
        self.mk_soma(name = 'myFirstSoma', diam = 5, L = 100)
        - self.soma becomes section with name = myFirstSoma, diam = 5 and L = 100
        """
        
        if name:
            self.soma = h.Section(name = name)
        else:
            self.soma = h.Section(name = 'soma')
        
        if diam:
            self.soma.diam = diam
        if L:
            self.soma = L = L
        
        return self.soma
    
    def mk_dendrite(self, n_secs, names = None, diam = None, L = None):
        """Adds a dendrite to list self.dendrites. Before first call
        self.dendrites = None. On first call self.dendrite becomes list with 1
        dendrite. Automatically connects the first section of the dendrite to
        the soma. Raises error if self.soma = None.
        """
        if not self.dendrites:
            self.dendrites = []
            
        
        
    def mk_sections(self, num_dend, num_seg):
        """Sets up the Sections WITH Topology"""
        self.soma = h.Section(name='soma')
        self.all_sections.append(self.soma)
        self.dendrites = []
        self.num_dend = num_dend
        self.i = 0  # Counter for __iter__

        if hasattr(num_seg, '__iter'):
            for curr_dend in range(num_dend):
                for curr_num_segs in num_seg:
                    self.dendrites.append(Dendrite(n_segs = curr_num_segs))
        else:
            for curr_dend in range(num_dend):
                self.dendrites.append(Dendrite(n_segs = num_seg))
        for x in self.dendrites:
            for y in x:
                self.all_sections.append(y)

        for x in self.dendrites:
            x.connect_segments(self.soma)

    def mk_geometry(self, soma_diam, soma_L, dend_L, dend_diam):
        """Sets up the geometry of the sections
        dend_L - scalar or iterable matching the number of dendrites
        dend_diam - scalar or iterable matching the number of dendrites
        """
        self.soma.diam = soma_diam
        self.soma.L = soma_L

        if len(dend_L) != self.num_dend:
            raise ValueError("dend_L must match the number of dendrites or be a scalar")
        if len(dend_diam) != self.num_dend:
            raise ValueError("dend_diam must match the number of dendrites or be a scalar")
        for dend_idx, x in enumerate(self.dendrites):
            if len(self.dendrites[dend_idx]) != len(dend_L[dend_idx]):
                raise ValueError
            if len(self.dendrites[dend_idx]) != len(dend_diam[dend_idx]):
                raise ValueError
            for seg_idx, y in enumerate(x):
                y.L = dend_L[dend_idx][seg_idx]
                y.diam = dend_diam[dend_idx][seg_idx]

    def connect_post(self, source, synapse, thr = 10, delay = 1, weight = 0):
        if not hasattr(self, 'synapses'):
            self.synapses = []
        if not hasattr(self, 'netcons'):
            self.netcons = []

        self.synapses.append(synapse)
        netcon = h.NetCon(source, synapse, thr, delay, weight)
        self.netcons.append(netcon)

        return netcon

    def _current_clamp_soma(self, amp = 0.3, dur = 500, delay = 500):
        """Setup a current clamp recording"""
        self.stim = h.IClamp(self.soma(1))
        self.stim.amp = amp     #Too high amps crash the simulation without raising an error!
        self.stim.dur = dur
        self.stim.delay = delay

        return self._voltage_recording()

    def _voltage_recording(self):
        soma_v_vec = h.Vector()
        t_vec = h.Vector()
        soma_v_vec.record(self.soma(1)._ref_v)
        t_vec.record(h._ref_t)
        
        return soma_v_vec, t_vec

    def __iter__(self):
        return self

    def next(self):
        if self.i < (len(self.segs)):
            i = self.i
            self.i += 1
            return self.all_sections[i]
        else:
            self.i = 0
            raise StopIteration()