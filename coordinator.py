# coding: utf-8

import math
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

        last_reactions = interrupter(self.sim).last_reactions()
        dirty = False
        for rr in last_reactions:
            if self._interrupt(t, rr[1]):
                dirty = True
        if dirty:
            self.sim.initialize()
            return True
        return False

class DiscreteEvent(SimulatorEvent):

    def __init__(self, sim):
        SimulatorEvent.__init__(self, sim)

    def next_time(self):
        return self.sim.next_time()

    def step(self):
        self.sim.step()

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

class ODESimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs, last_reactions=[], rng=GSLRandomNumberGenerator()):
        SimulatorAdapter.__init__(self, lhs, rhs)
        self.rng = rng
        self.__last_reactions = last_reactions
        assert isinstance(self.lhs, ode.ODESimulator)

    def last_reactions(self):
        if isinstance(self.rhs, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return self.__last_reactions
        elif isinstance(self.rhs, meso.MesoscopicSimulator):
            def convert(ri):
                return meso.ReactionInfo(
                    ri.t(), ri.reactants(), ri.products(),
                    self.rng.uniform_int(
                        0, self.rhs.world().num_subvolumes() - 1))
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        elif isinstance(self.rhs, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                coord = self.rng.uniform_int(0, self.rhs.world().size() - 1)
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        elif isinstance(self.rhs, egfrd.EGFRDSimulator):
            def convert(ri):
                lengths = self.lhs.world().edge_lengths()
                pos = Real3(self.rng.uniform(0, 1) * lengths[0],
                            self.rng.uniform(0, 1) * lengths[1],
                            self.rng.uniform(0, 1) * lengths[2])
                reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        raise ValueError("Not supported yet.")

class ODEEvent(DiscreteTimeEvent):

    def __init__(self, sim, dt=1.0):
        DiscreteTimeEvent.__init__(self, sim, dt)
        self.__last_reactions = []

    def step(self):
        DiscreteTimeEvent.step(self)
        self.__last_reactions = self.generate_reactions()

    def sync(self):
        self.__last_reactions = []

        dirty = False
        for sp in self.sim.world().list_species():
            if own(self.sim, sp):
                continue
            value = self.sim.world().get_value_exact(sp)
            if value >= 1.0:
                self.sim.world().set_value(sp, value - math.floor(value))
                dirty = True

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [sp for sp in ri.products() if own(self.sim, sp)]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.sim.world().add_molecules(sp, 1)
        return True

    def __call__(self, rhs):
        assert self.sim != rhs
        return ODESimulatorAdapter(self.sim, rhs, self.__last_reactions)

    def generate_reactions(self):
        retval = []
        for sp in self.sim.world().list_species():
            if own(self.sim, sp):
                continue
            value = self.sim.world().get_value_exact(sp)
            if value >= 1.0:
                rr = ReactionRule([], [sp], 0.0)
                ri = gillespie.ReactionInfo(self.sim.t(), [], [sp])
                for i in range(int(value)):
                    retval.append((rr, ri))
        return retval

class GillespieWorldAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def remove_molecules(self, sp, num, coord=None):
        self.lhs.remove_molecules(sp, num)

    def remove_voxel(self, voxel):
        pid, v = voxel
        self.lhs.remove_molecules(v.species(), 1)

    def __getattr__(self, name):
        return getattr(self.lhs, name)

class GillespieSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs, gillespie.GillespieSimulator)

    def last_reactions(self):
        if isinstance(self.rhs, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return self.lhs.last_reactions()
        elif isinstance(self.rhs, meso.MesoscopicSimulator):
            def convert(ri):
                return meso.ReactionInfo(
                    ri.t(), ri.reactants(), ri.products(),
                    self.lhs.world().rng().uniform_int(
                        0, self.rhs.world().num_subvolumes() - 1))
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                coord = self.lhs.world().rng().uniform_int(0, self.rhs.world().size() - 1)
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, egfrd.EGFRDSimulator):
            def convert(ri):
                rng = self.lhs.world().rng()
                lengths = self.lhs.world().edge_lengths()
                pos = Real3(rng.uniform(0, 1) * lengths[0],
                            rng.uniform(0, 1) * lengths[1],
                            rng.uniform(0, 1) * lengths[2])
                reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        raise ValueError("Not supported yet.")

    def world(self):
        return GillespieWorldAdapter(self.lhs.world(), self.rhs.world())

class GillespieEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]

        dirty = False
        for sp in ri.products():
            if own(self.sim, sp):
                continue
            dirty = True
            self.sim.world().remove_molecules(sp, 1)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [sp for sp in ri.products() if own(self.sim, sp)]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.sim.world().add_molecules(sp, 1)
        return True

    def __call__(self, rhs):
        assert self.sim != rhs
        return GillespieSimulatorAdapter(self.sim, rhs)

class MesoscopicWorldAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    def remove_voxel(self, voxel):
        pid, v = voxel
        pos = self.rhs.coordinate2position(v.coordinate())
        coord = self.lhs.position2coordinate(pos)
        assert self.lhs.num_molecules(v.species(), coord) > 0
        self.lhs.remove_molecules(v.species(), 1, coord)

    def __getattr__(self, name):
        return getattr(self.lhs, name)

class MesoscopicSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs, meso.MesoscopicSimulator)

    def last_reactions(self):
        if isinstance(self.rhs, meso.MesoscopicSimulator):
            return self.lhs.last_reactions()
        elif isinstance(self.rhs, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return [(rr, gillespie.ReactionInfo(ri.t(), ri.reactants(), ri.products()))
                    for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                g = self.lhs.world().coord2global(ri.coordinate())
                rng = self.lhs.world().rng()
                lengths = self.lhs.world().subvolume_edge_lengths()
                pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                            (g[1] + rng.uniform(0, 1)) * lengths[1],
                            (g[2] + rng.uniform(0, 1)) * lengths[2])
                assert ri.coordinate() == self.lhs.world().position2coordinate(pos)
                coord = self.rhs.world().position2coordinate(pos)  #XXX: This may cause an overlap
                # assert ri.coordinate() == self.lhs.world().position2coordinate(self.rhs.world().coordinate2position(coord))  #XXX: This is not always True.
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, egfrd.EGFRDSimulator):
            def convert(ri):
                g = self.lhs.world().coord2global(ri.coordinate())
                rng = self.lhs.world().rng()
                lengths = self.lhs.world().subvolume_edge_lengths()
                pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                            (g[1] + rng.uniform(0, 1)) * lengths[1],
                            (g[2] + rng.uniform(0, 1)) * lengths[2])
                reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs)))

    # def world(self):
    #     return MesoscopicWorldAdapter(self.lhs.world(), self.rhs.world())

class MesoscopicEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]
        coord = ri.coordinate()

        dirty = False
        for sp in ri.products():
            if own(self.sim, sp):
                continue
            dirty = True
            self.sim.world().remove_molecules(sp, 1, coord)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [sp for sp in ri.products() if own(self.sim, sp)]
        if len(products) == 0:
            return False
        coord = ri.coordinate()
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.sim.world().add_molecules(sp, 1, coord)
        return True

    def __call__(self, rhs):
        assert self.sim != rhs
        return MesoscopicSimulatorAdapter(self.sim, rhs)

class SpatiocyteWorldAdapter:

    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs

    #def remove_molecules(self, sp, num, coord=None):
    #    if coord is None:
    #        self.lhs.remove_molecules(sp, num)
    #        return
    #    pids = [pid for pid, v in self.lhs.list_voxels_exact(sp)
    #            if self.rhs.position2coordinate(
    #                self.lhs.coordinate2position(v.coordinate())) == coord]
    #    assert len(pids) >= num
    #    for i in range(num):
    #        pid = pids.pop(self.lhs.rng().uniform_int(0, len(pids) - 1))
    #        self.lhs.remove_voxel(pid)

    def list_coordinates_exact(self, sp):
        return [self.rhs.position2coordinate(
                    self.lhs.coordinate2position(v.coordinate()))
                for pid, v in self.lhs.list_voxels_exact(sp)]

    def __getattr__(self, name):
        return getattr(self.lhs, name)

class SpatiocyteSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs, spatiocyte.SpatiocyteSimulator)

    def last_reactions(self):
        if isinstance(self.rhs, spatiocyte.SpatiocyteSimulator):
            return self.lhs.last_reactions()
        elif isinstance(self.rhs, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = [v.species() for pid, v in ri.reactants()]
                products = [v.species() for pid, v in ri.products()]
                return gillespie.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = [v.species() for pid, v in ri.reactants()]
                products = [v.species() for pid, v in ri.products()]
                assert len(products) > 0
                pos = self.lhs.world().coordinate2position(ri.products()[0][1].coordinate())
                coord = self.rhs.world().position2coordinate(pos)
                return meso.ReactionInfo(ri.t(), reactants, products, coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, egfrd.EGFRDSimulator):
            def convert(ri):
                reactants = [(pid, Particle(v.species(), self.lhs.world().coordinate2position(v.coordinate()), v.radius(), v.D()))
                              for pid, v in ri.reactants()]
                products = [(pid, Particle(v.species(), self.lhs.world().coordinate2position(v.coordinate()), v.radius(), v.D()))
                             for pid, v in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs)))

    # def world(self):
    #     return SpatiocyteWorldAdapter(self.lhs.world(), self.rhs.world())

class SpatiocyteEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]

        dirty = False
        for pid, v in ri.products():
            if own(self.sim, v.species()):
                continue
            dirty = True
            self.sim.world().remove_voxel(pid)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [(pid, v) for pid, v in ri.products() if own(self.sim, v.species())]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for pid, v in products:
            self.sim.world().new_voxel(v.species(), v.coordinate())
        return True

    def __call__(self, rhs):
        assert self.sim != rhs
        return SpatiocyteSimulatorAdapter(self.sim, rhs)

class EGFRDSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs, egfrd.EGFRDSimulator)

    def last_reactions(self):
        if isinstance(self.rhs, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                reactants = [(pid, Voxel(v.species(), self.rhs.world().position2coordinate(v.position()), v.radius(), v.D()))
                              for pid, v in ri.reactants()]
                products = [(pid, Voxel(v.species(), self.rhs.world().position2coordinate(v.position()), v.radius(), v.D()))
                             for pid, v in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = [p.species() for pid, p in ri.reactants()]
                products = [p.species() for pid, p in ri.products()]
                return gillespie.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = [p.species() for pid, p in ri.reactants()]
                products = [p.species() for pid, p in ri.products()]
                assert len(products) > 0
                coord = self.rhs.world().position2coordinate(ri.products()[0][1].position())
                return meso.ReactionInfo(ri.t(), reactants, products, coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.last_reactions()]
        elif isinstance(self.rhs, egfrd.EGFRDSimulator):
            return self.lhs.last_reactions()
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs)))

class EGFRDEvent(DiscreteEvent):

    def __init__(self, sim):
        DiscreteEvent.__init__(self, sim)

    def sync(self):
        last_reactions = self.sim.last_reactions()
        if len(last_reactions) == 0:
            return
        assert len(last_reactions) == 1
        ri = last_reactions[0][1]

        dirty = False
        for pid, p in ri.products():
            if own(self.sim, p.species()):
                continue
            dirty = True
            self.sim.world().remove_particle(pid)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [(pid, p) for pid, p in ri.products() if own(self.sim, p.species())]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for pid, p in products:
            self.sim.world().new_particle(p.species(), p.position())
        return True

    def __call__(self, rhs):
        assert self.sim != rhs
        return EGFRDSimulatorAdapter(self.sim, rhs)

def own(sim, sp):
    if sp in (Species("A1"), Species("A2")):
        return isinstance(sim, gillespie.GillespieSimulator)
    elif sp in (Species("B1"), Species("B2")):
        return isinstance(sim, meso.MesoscopicSimulator)
    elif sp in (Species("C1"), Species("C2")):
        return isinstance(sim, spatiocyte.SpatiocyteSimulator)
    elif sp in (Species("D1"), Species("D2")):
        return isinstance(sim, egfrd.EGFRDSimulator)
    elif sp in (Species("E1"), Species("E2")):
        return isinstance(sim, ode.ODESimulator)
    raise ValueError("Unknown species [{}] was given".format(sp.serial()))

class Coordinator:

    def __init__(self):
        self.events = []
        self.num_steps = 0
        self.last_event = None

    def add_simulator(self, sim):
        if isinstance(sim, ode.ODESimulator):
            ev = ODEEvent(sim)
        elif isinstance(sim, gillespie.GillespieSimulator):
            ev = GillespieEvent(sim)
        elif isinstance(sim, meso.MesoscopicSimulator):
            ev = MesoscopicEvent(sim)
        elif isinstance(sim, spatiocyte.SpatiocyteSimulator):
            ev = SpatiocyteEvent(sim)
        elif isinstance(sim, egfrd.EGFRDSimulator):
            ev = EGFRDEvent(sim)
        else:
            raise ValueError("{} not supported.".format(repr(sim)))

        self.add_event(ev)

    def add_event(self, ev):
        self.events.append(ev)

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

        self.interrupt_all(ntime, self.last_event, (idx, ))

        self.last_event.sync()  #XXX: under development
