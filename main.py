# coding: utf-8

import numpy
from ecell4 import *


D, radius = 1, 0.005
edge_lengths = Real3(1, 1, 1)

with species_attributes():
    A1 | A2 | B1 | B2 | C1 | C2 | D1 | D2 | {"D": str(D), "radius": str(radius)}

with reaction_rules():
    A1 == A2 | (1.0, 1.0)
    B1 == B2 | (1.0, 1.0)
    C1 == C2 | (1.0, 1.0)
    D1 == D2 | (1.0, 1.0)

    A1 == B1 | (1.0, 1.0)
    B1 == C1 | (1.0, 1.0)
    C1 == A1 | (1.0, 1.0)
    A1 == D1 | (1.0, 1.0)
    B1 == D1 | (1.0, 1.0)
    C1 == D1 | (1.0, 1.0)

m = get_model()

w1 = gillespie.GillespieWorld(edge_lengths)
w1.bind_to(m)
sim1 = gillespie.GillespieSimulator(w1)

w2 = meso.MesoscopicWorld(edge_lengths, Integer3(9, 9, 9))
w2.bind_to(m)
sim2 = meso.MesoscopicSimulator(w2)

w3 = spatiocyte.SpatiocyteWorld(edge_lengths, radius)
w3.bind_to(m)
sim3 = spatiocyte.SpatiocyteSimulator(w3)

w4 = egfrd.EGFRDWorld(edge_lengths, Integer3(4, 4, 4))
w4.bind_to(m)
sim4 = egfrd.EGFRDSimulator(m, w4)

w1.add_molecules(Species("A1"), 240)
# w1.add_molecules(Species("A1"), 60)
# w2.add_molecules(Species("B1"), 60)
# w3.add_molecules(Species("C1"), 60)
# w4.add_molecules(Species("D1"), 60)

import coordinator
owner = coordinator.Coordinator()
owner.add_simulator(sim1)
owner.add_simulator(sim2)
owner.add_simulator(sim3)
owner.add_simulator(sim4)
owner.initialize()

data = []
def log(owner):
    data.append((
        owner.last_event.t(),
        owner.events[0].sim.world().num_molecules_exact(Species("A1")),
        owner.events[0].sim.world().num_molecules_exact(Species("A2")),
        owner.events[1].sim.world().num_molecules_exact(Species("B1")),
        owner.events[1].sim.world().num_molecules_exact(Species("B2")),
        owner.events[2].sim.world().num_molecules_exact(Species("C1")),
        owner.events[2].sim.world().num_molecules_exact(Species("C2")),
        owner.events[3].sim.world().num_molecules_exact(Species("D1")),
        owner.events[3].sim.world().num_molecules_exact(Species("D2")),
        ))

while owner.step(3):
    if len(owner.last_event.sim.last_reactions()) > 0:
        log(owner)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

data = numpy.asarray(data)
for i in range(1, 9):
    plt.plot(data.T[0], data.T[i], '-')
# plt.show()
plt.savefig('result.eps')

for ev in owner.events:
    print('=> {}, {}'.format(ev.t(), ev.num_steps()))
    print([ev.sim.world().num_molecules_exact(Species(name))
           for name in ("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2", "C1_")])
