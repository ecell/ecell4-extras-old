# coding: utf-8

import math
import collections
from coordinator import SimulatorAdapter, DiscreteEvent, DiscreteTimeEvent
# from ecell4 import *
from ecell4.core import Real3, Integer3, Species, ReactionRule, GSLRandomNumberGenerator, ParticleID, Voxel, Particle
from ecell4 import gillespie, meso, spatiocyte, egfrd, ode


def simulator_event(sim):
    if isinstance(sim, ode.ODESimulator):
        return ODEEvent(sim)
    elif isinstance(sim, gillespie.GillespieSimulator):
        return GillespieEvent(sim)
    elif isinstance(sim, meso.MesoscopicSimulator):
        return MesoscopicEvent(sim)
    elif isinstance(sim, spatiocyte.SpatiocyteSimulator):
        return SpatiocyteEvent(sim)
    elif isinstance(sim, egfrd.EGFRDSimulator):
        return EGFRDEvent(sim)
    raise ValueError("{} not supported.".format(repr(sim)))

class ODESimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs, last_reactions=[], rng=GSLRandomNumberGenerator()):
        SimulatorAdapter.__init__(self, lhs, rhs)
        self.rng = rng
        self.__last_reactions = last_reactions
        assert isinstance(self.lhs.sim, ode.ODESimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return self.__last_reactions
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                return meso.ReactionInfo(
                    ri.t(), ri.reactants(), ri.products(),
                    self.rng.uniform_int(
                        0, self.rhs.sim.world().num_subvolumes() - 1))
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        elif isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                coord = self.rng.uniform_int(0, self.rhs.sim.world().size() - 1)
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.__last_reactions]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                lengths = self.lhs.sim.world().edge_lengths()
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

    def __init__(self, sim, dt=None):
        if dt is None:
            dt = sim.dt()
        DiscreteTimeEvent.__init__(self, sim, dt)
        self.__last_reactions = []

    def step(self):
        DiscreteTimeEvent.step(self)
        self.__last_reactions = self.generate_reactions()

    def sync(self):
        self.__last_reactions = []

        dirty = False
        for sp in self.sim.world().list_species():
            if self.own(sp):
                continue
            value = self.sim.world().get_value_exact(sp)
            if value >= 1.0:
                self.sim.world().set_value(sp, value - math.floor(value))
                dirty = True

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [sp for sp in ri.products() if self.own(sp)]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.sim.world().add_molecules(sp, 1)
        return True

    def __call__(self, rhs):
        assert self.sim != rhs.sim
        return ODESimulatorAdapter(self, rhs, self.__last_reactions)

    def generate_reactions(self):
        retval = []
        for sp in self.sim.world().list_species():
            if self.own(sp):
                continue
            value = self.sim.world().get_value_exact(sp)
            if value >= 1.0:
                rr = ReactionRule([], [sp], 0.0)
                ri = gillespie.ReactionInfo(self.sim.t(), [], [sp])
                for i in range(int(value)):
                    retval.append((rr, ri))
        return retval

    def updated(self):
        return len(self.__last_reactions) > 0

# class GillespieWorldAdapter:
# 
#     def __init__(self, lhs, rhs):
#         self.lhs = lhs
#         self.rhs = rhs
# 
#     def remove_molecules(self, sp, num, coord=None):
#         self.lhs.remove_molecules(sp, num)
# 
#     def remove_voxel(self, voxel):
#         pid, v = voxel
#         self.lhs.remove_molecules(v.species(), 1)
# 
#     def __getattr__(self, name):
#         return getattr(self.lhs, name)

class GillespieSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, gillespie.GillespieSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return self.lhs.sim.last_reactions()
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                return meso.ReactionInfo(
                    ri.t(), ri.reactants(), ri.products(),
                    self.lhs.sim.world().rng().uniform_int(
                        0, self.rhs.sim.world().num_subvolumes() - 1))
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                coord = self.lhs.sim.world().rng().uniform_int(0, self.rhs.sim.world().size() - 1)
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                rng = self.lhs.sim.world().rng()
                lengths = self.lhs.sim.world().edge_lengths()
                pos = Real3(rng.uniform(0, 1) * lengths[0],
                            rng.uniform(0, 1) * lengths[1],
                            rng.uniform(0, 1) * lengths[2])
                reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        raise ValueError("Not supported yet.")

    # def world(self):
    #     return GillespieWorldAdapter(self.lhs.world(), self.rhs.world())

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
            if self.own(sp):
                continue
            dirty = True
            self.sim.world().remove_molecules(sp, 1)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [sp for sp in ri.products() if self.own(sp)]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.sim.world().add_molecules(sp, 1)
        return True

    def _mirror(self, interrupter, src, dst):
        w = interrupter.sim.world()
        value1 = w.get_value_exact(src)
        value2 = self.sim.world().get_value_exact(dst)
        if value1 != value2:
            self.sim.world().set_value(dst, value1)
            return True
        return False

    def __call__(self, rhs):
        assert self.sim != rhs.sim
        return GillespieSimulatorAdapter(self, rhs)

# class MesoscopicWorldAdapter:
# 
#     def __init__(self, lhs, rhs):
#         self.lhs = lhs
#         self.rhs = rhs
# 
#     def remove_voxel(self, voxel):
#         pid, v = voxel
#         pos = self.rhs.coordinate2position(v.coordinate())
#         coord = self.lhs.position2coordinate(pos)
#         assert self.lhs.num_molecules(v.species(), coord) > 0
#         self.lhs.remove_molecules(v.species(), 1, coord)
# 
#     def __getattr__(self, name):
#         return getattr(self.lhs, name)

class MesoscopicSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, meso.MesoscopicSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            return self.lhs.sim.last_reactions()
        elif isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            return [(rr, gillespie.ReactionInfo(ri.t(), ri.reactants(), ri.products()))
                    for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                g = self.lhs.sim.world().coord2global(ri.coordinate())
                rng = self.lhs.sim.world().rng()
                lengths = self.lhs.sim.world().subvolume_edge_lengths()
                pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                            (g[1] + rng.uniform(0, 1)) * lengths[1],
                            (g[2] + rng.uniform(0, 1)) * lengths[2])
                assert ri.coordinate() == self.lhs.sim.world().position2coordinate(pos)
                coord = self.rhs.sim.world().position2coordinate(pos)  #XXX: This may cause an overlap
                # assert ri.coordinate() == self.lhs.sim.world().position2coordinate(self.rhs.sim.world().coordinate2position(coord))  #XXX: This is not always True.
                reactants = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Voxel(sp, coord, 0, 0))
                             for sp in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                g = self.lhs.sim.world().coord2global(ri.coordinate())
                rng = self.lhs.sim.world().rng()
                lengths = self.lhs.sim.world().subvolume_edge_lengths()
                pos = Real3((g[0] + rng.uniform(0, 1)) * lengths[0],
                            (g[1] + rng.uniform(0, 1)) * lengths[1],
                            (g[2] + rng.uniform(0, 1)) * lengths[2])
                reactants = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.reactants()]
                products = [(ParticleID(), Particle(sp, pos, 0, 0))
                             for sp in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs.sim)))

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
            if self.own(sp):
                continue
            dirty = True
            self.sim.world().remove_molecules(sp, 1, coord)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [sp for sp in ri.products() if self.own(sp)]
        if len(products) == 0:
            return False
        coord = ri.coordinate()
        self.sim.step(t)
        assert self.sim.t() == t
        for sp in products:
            self.sim.world().add_molecules(sp, 1, coord)
        return True

    def __call__(self, rhs):
        assert self.sim != rhs.sim
        return MesoscopicSimulatorAdapter(self, rhs)

# class SpatiocyteWorldAdapter:
# 
#     def __init__(self, lhs, rhs):
#         self.lhs = lhs
#         self.rhs = rhs
# 
#     #def remove_molecules(self, sp, num, coord=None):
#     #    if coord is None:
#     #        self.lhs.remove_molecules(sp, num)
#     #        return
#     #    pids = [pid for pid, v in self.lhs.list_voxels_exact(sp)
#     #            if self.rhs.position2coordinate(
#     #                self.lhs.coordinate2position(v.coordinate())) == coord]
#     #    assert len(pids) >= num
#     #    for i in range(num):
#     #        pid = pids.pop(self.lhs.rng().uniform_int(0, len(pids) - 1))
#     #        self.lhs.remove_voxel(pid)
# 
#     def list_coordinates_exact(self, sp):
#         return [self.rhs.position2coordinate(
#                     self.lhs.coordinate2position(v.coordinate()))
#                 for pid, v in self.lhs.list_voxels_exact(sp)]
# 
#     def __getattr__(self, name):
#         return getattr(self.lhs, name)

class SpatiocyteSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, spatiocyte.SpatiocyteSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            return self.lhs.sim.last_reactions()
        elif isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = [v.species() for pid, v in ri.reactants()]
                products = [v.species() for pid, v in ri.products()]
                return gillespie.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = [v.species() for pid, v in ri.reactants()]
                products = [v.species() for pid, v in ri.products()]
                assert len(products) > 0
                pos = self.lhs.sim.world().coordinate2position(ri.products()[0][1].coordinate())
                coord = self.rhs.sim.world().position2coordinate(pos)
                return meso.ReactionInfo(ri.t(), reactants, products, coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            def convert(ri):
                reactants = [(pid, Particle(v.species(), self.lhs.sim.world().coordinate2position(v.coordinate()), v.radius(), v.D()))
                              for pid, v in ri.reactants()]
                products = [(pid, Particle(v.species(), self.lhs.sim.world().coordinate2position(v.coordinate()), v.radius(), v.D()))
                             for pid, v in ri.products()]
                return egfrd.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs.sim)))

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
            if self.own(v.species()):
                continue
            dirty = True
            self.sim.world().remove_voxel(pid)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [(pid, v) for pid, v in ri.products() if self.own(v.species())]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for pid, v in products:
            self.sim.world().new_voxel(v.species(), v.coordinate())
        return True

    def __call__(self, rhs):
        assert self.sim != rhs.sim
        return SpatiocyteSimulatorAdapter(self, rhs)

class EGFRDSimulatorAdapter(SimulatorAdapter):

    def __init__(self, lhs, rhs):
        SimulatorAdapter.__init__(self, lhs, rhs)
        assert isinstance(self.lhs.sim, egfrd.EGFRDSimulator)

    def last_reactions(self):
        if isinstance(self.rhs.sim, spatiocyte.SpatiocyteSimulator):
            def convert(ri):
                reactants = [(pid, Voxel(v.species(), self.rhs.sim.world().position2coordinate(v.position()), v.radius(), v.D()))
                              for pid, v in ri.reactants()]
                products = [(pid, Voxel(v.species(), self.rhs.sim.world().position2coordinate(v.position()), v.radius(), v.D()))
                             for pid, v in ri.products()]
                return spatiocyte.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, (ode.ODESimulator, gillespie.GillespieSimulator)):
            def convert(ri):
                reactants = [p.species() for pid, p in ri.reactants()]
                products = [p.species() for pid, p in ri.products()]
                return gillespie.ReactionInfo(ri.t(), reactants, products)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, meso.MesoscopicSimulator):
            def convert(ri):
                reactants = [p.species() for pid, p in ri.reactants()]
                products = [p.species() for pid, p in ri.products()]
                assert len(products) > 0
                coord = self.rhs.sim.world().position2coordinate(ri.products()[0][1].position())
                return meso.ReactionInfo(ri.t(), reactants, products, coord)
            return [(rr, convert(ri)) for (rr, ri) in self.lhs.sim.last_reactions()]
        elif isinstance(self.rhs.sim, egfrd.EGFRDSimulator):
            return self.lhs.sim.last_reactions()
        raise ValueError("Not supported yet [{}].".format(repr(self.rhs.sim)))

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
            if self.own(p.species()):
                continue
            dirty = True
            self.sim.world().remove_particle(pid)

        if dirty:
            self.sim.initialize()

    def _interrupt(self, t, ri):
        products = [(pid, p) for pid, p in ri.products() if self.own(p.species())]
        if len(products) == 0:
            return False
        self.sim.step(t)
        assert self.sim.t() == t
        for pid, p in products:
            self.sim.world().new_particle(p.species(), p.position())
        return True

    def __call__(self, rhs):
        assert self.sim != rhs.sim
        return EGFRDSimulatorAdapter(self, rhs)
