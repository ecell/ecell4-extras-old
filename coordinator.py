# coding: utf-8

import math
import collections
from ecell4 import *


class SimulatorAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def __getattr__(self, name):
        return getattr(self.lhs, name)

class SimulatorEvent:

    def __init__(self, sim):
        self.sim = sim
        self.__belongings = set()
        self.__borrowings = dict()

    def add(self, value):
        if isinstance(value, Species):
            self.__belongings.add(value)
        elif isinstance(value, str):
            self.__belongings.add(Species(value))
        elif isinstance(value, collections.Iterable):
            for sp in value:
                self.add(sp)
        else:
            raise ValueError(
                'an invalid argument [{}] was given.'.format(repr(value)))

    def borrow(self, src, dst):
        if isinstance(src, Species) and isinstance(dst, Species):
            self.__borrowings[dst] = src
        elif isinstance(src, str) and isinstance(dst, str):
            self.__borrowings[Species(dst)] = Species(src)
        else:
            raise ValueError(
                'an invalid argument [{} and {}] was given.'.format(
                    repr(src), repr(dst)))

    def own(self, sp):
        return (sp in self.__belongings)

    def owe(self, sp):
        return self.__borrowings.get(sp, None)

    def initialize(self):
        self.sim.initialize()

    def t(self):
        return self.sim.t()

    def num_steps(self):
        return self.sim.num_steps()

    def interrupt(self, t, interrupter):
        if interrupter is None:
            self.sim.step(t)
            assert self.sim.t() == t
            return True

        dirty = False

        for dst, src in self.__borrowings.items():
            if not interrupter.own(src):
                continue
            if self._mirror(interrupter, src, dst):
                dirty = True

        if interrupter.updated():
            last_reactions = interrupter(self.sim).last_reactions()
            for rr in last_reactions:
                if self._interrupt(t, rr[1]):
                    dirty = True

        if dirty:
            self.sim.initialize()
        return dirty

    def updated(self):
        return True

    def _mirror(self, interrupter, src, dst):
        raise RuntimeError('Not implemented yet [{}].'.format(repr(self)))

class DiscreteEvent(SimulatorEvent):

    def __init__(self, sim):
        SimulatorEvent.__init__(self, sim)

    def next_time(self):
        return self.sim.next_time()

    def step(self):
        self.sim.step()

    def updated(self):
        return len(self.sim.last_reactions()) > 0

class DiscreteTimeEvent(SimulatorEvent):

    def __init__(self, sim, dt):
        SimulatorEvent.__init__(self, sim)
        self.__t = 0.0
        self.__dt = dt
        self.__num_steps = 0

    def next_time(self):
        return self.__t + self.__dt * (self.__num_steps + 1)

    def step(self):
        assert self.sim.t() <= self.next_time()
        self.sim.step(self.next_time())
        self.__num_steps += 1

class Coordinator:

    def __init__(self):
        self.events = []
        self.__t = 0.0
        self.num_steps = 0
        self.last_event = None

    def t(self):
        return self.__t

    def get_value(self, sp):
        value = 0.0
        for ev in self.events:
            if ev.own(sp):
                value += ev.sim.world().get_value_exact(sp)
        return value

    def add_event(self, ev):
        self.events.append(ev)
        return ev

    def initialize(self):
        for ev in self.events:
            ev.initialize()
        self.last_event = None

    def get_next_event(self):
        idx = 0
        ntime = self.events[idx].next_time()
        for i, ev in enumerate(self.events[1:]):
            if ev.next_time() < ntime:
                (idx, ntime) = (i + 1, ev.next_time())
        return (idx, ntime)

    def interrupt_all(self, t, interrupter=None, ignore=()):
        for i, ev in enumerate(self.events):
            if i in ignore:
                continue
            if ev.interrupt(t, interrupter):
                self.interrupt_all(t, ev, (i, ))

    def step(self, upto=None):
        self.num_steps += 1
        self.last_event = None

        if upto is None:
            self._step()
            return False

        idx, ntime = self.get_next_event()
        if ntime > upto:
            self.__t = upto
            self.interrupt_all(upto, None)
            return False

        self._step()
        return True

    def _step(self):
        if len(self.events) == 0:
            return

        idx, ntime = self.get_next_event()
        self.last_event = self.events[idx]
        self.last_event.step()
        assert self.last_event.t() == ntime
        self.__t = ntime

        self.interrupt_all(ntime, self.last_event, (idx, ))

        self.last_event.sync()  #XXX: under development
