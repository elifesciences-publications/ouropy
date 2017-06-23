# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 12:44:23 2017

@author: DanielM
"""
from neuron import h, gui

class GenDendrite(object):
    """This is the model of a generic dendrite."""
    def __init__(self, dend_name = None, sec_names = None, n_secs= None, diam = None, L = None):
        self.secs = None
        """if name:
            self.name = dend_name
        else:
            self.name = 'dendrite'
        if n_secs:
            if not (type(sec_names) == list):
                
        self.secs = None
        self.i = 0 #initializing i for __iter__ protocol
        self.mk_segments(n_secs = n_secs)
        if bool(diam):
            self.set_diam(diam)
        if bool(L):
            self.set_L(L)"""

    def mk_secs(self, n_secs = 1, sec_names = None):
        """Makes sections AND connects them. This is because a dendrite is by
        definition made up of connected sections. sec_names has to be a list
        of section names with len = n_secs. If sec_names = None the section
        names are 'sec' + str(number).
        """
        if not self.secs:
            self.secs = []
        if sec_names:
            if not (type(sec_names) == list):
                raise TypeError("sec_names must be list or None")
            if len(sec_names) != n_secs:
                raise ValueError("The len of sec_names must equal n_secs")

        for curr_n in range(n_secs):
            if sec_names:
                self.secs.append(h.Section(name = sec_names[curr_n]))
            else:
                self.secs.append(h.Section(name = 'sec' + str(curr_n)))
            if curr_n > 0:
                self.secs[curr_n].connect(self.secs[curr_n - 1](1))

    def conn_soma(self):
        pass
    def set_diam(self, diam):
        if not bool(self.secs):
            raise Warning("Can't set diameter before segments are made")
            return
        if hasattr(diam, '__iter__'):
            if len(diam) != len(self.secs):
                raise Warning("List of diameters does not fit number of segments")
                return
            for idx, curr_seg in enumerate(self.secs):
                curr_seg.diam = diam[idx]
            return
        else:
            for curr_seg in self.secs:
                curr_seg.diam = diam

    def set_L(self, L):
        if not bool(self.secs):
            raise Warning("Can't set L before segments are made")
            return
        if hasattr(L, '__iter__'):
            if len(L) != len(self.secs):
                raise Warning("List of diameters does not fit number of segments")
                return
            for idx, curr_seg in enumerate(self.secs):
                curr_seg.L = L[idx]
            return
        else:
            for curr_seg in self.secs:
                curr_seg.L = L

    def __iter__(self):
        return self

    def next(self):
        if self.i < (len(self.secs)):
            i = self.i
            self.i += 1
            return self.secs[i]
        else:
            self.i = 0
            raise StopIteration()

    def __getitem__(self, key):
        return self.secs[key]

    def __len__(self):
        return len(self.secs)