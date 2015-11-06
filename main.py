# coding: utf-8

import numpy

from ecell4.core import Real3, Integer3, Species
from ecell4 import gillespie, meso, spatiocyte, egfrd, ode
from ecell4.util import species_attributes, reaction_rules, get_model

import coordinator
from implementations import simulator_event, ODEEvent


def test1():
    D, radius = 1, 0.005
    edge_lengths = Real3(1, 1, 1)

    with species_attributes():
        A1 | A2 | B1 | B2 | C1 | C2 | D1 | D2 | E1 | E2 | {
            "D": str(D), "radius": str(radius)}

    with reaction_rules():
        A1 == A2 | (1.0, 1.0)
        B1 == B2 | (1.0, 1.0)
        C1 == C2 | (1.0, 1.0)
        D1 == D2 | (1.0, 1.0)
        E1 == E2 | (1.0, 1.0)

        A1 == B1 | (1.0, 1.0)
        B1 == C1 | (1.0, 1.0)
        C1 == A1 | (1.0, 1.0)
        A1 == D1 | (1.0, 1.0)
        B1 == D1 | (1.0, 1.0)
        C1 == D1 | (1.0, 1.0)

        A1 == E1 | (1.0, 1.0)
        B1 == E1 | (1.0, 1.0)
        C1 == E1 | (1.0, 1.0)
        D1 == E1 | (1.0, 1.0)


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

    w5 = ode.ODEWorld(edge_lengths)
    w5.bind_to(m)
    sim5 = ode.ODESimulator(w5)

    w1.add_molecules(Species("A1"), 300)
    # w1.add_molecules(Species("A1"), 60)
    # w2.add_molecules(Species("B1"), 60)
    # w3.add_molecules(Species("C1"), 60)
    # w4.add_molecules(Species("D1"), 60)
    # w5.add_molecules(Species("E1"), 60)

    owner = coordinator.Coordinator()
    owner.add_event(simulator_event(sim1)).register_species(('A1', 'A2'))
    owner.add_event(simulator_event(sim2)).register_species(('B1', 'B2'))
    owner.add_event(simulator_event(sim3)).register_species(('C1', 'C2'))
    owner.add_event(simulator_event(sim4)).register_species(('D1', 'D2'))
    owner.add_event(ODEEvent(sim5, 0.01)).register_species(('E1', 'E2'))
    owner.initialize()

    data = []
    def log(owner):
        data.append((
            owner.t(),
            owner.get_value(Species("A1")),
            owner.get_value(Species("A2")),
            owner.get_value(Species("B1")),
            owner.get_value(Species("B2")),
            owner.get_value(Species("C1")),
            owner.get_value(Species("C2")),
            owner.get_value(Species("D1")),
            owner.get_value(Species("D2")),
            owner.get_value(Species("E1")),
            owner.get_value(Species("E2")),
            ))

    log(owner)
    while owner.step(3):
        if owner.last_event.updated():
            log(owner)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt

    data = numpy.asarray(data)
    for i in range(1, 11):
        plt.plot(data.T[0], data.T[i], '-')
    # plt.show()
    plt.savefig('result.eps')

    # for ev in owner.events:
    #     print('=> {}, {}'.format(ev.t(), ev.num_steps()))
    #     print([ev.sim.world().num_molecules_exact(Species(name))
    #            for name in ("A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2", "C1_")])
    numpy.savetxt("result.txt", data)

def test2():
    edge_lengths = Real3(1, 1, 1)

    with reaction_rules():
        A1 == A2 | (1.0, 1.0)
        E1 == E2 | (1.0, 1.0)
        A1 > E1 | 1.0
        E1 > A1 | 1.0

    m = get_model()

    w1 = gillespie.GillespieWorld(edge_lengths)
    w1.bind_to(m)
    sim1 = gillespie.GillespieSimulator(w1)

    w2 = ode.ODEWorld(edge_lengths)
    w2.bind_to(m)
    sim2 = ode.ODESimulator(w2)

    w1.add_molecules(Species("A1"), 120)
    # w1.add_molecules(Species("A1"), 60)
    # w2.add_molecules(Species("E1"), 60)

    owner = coordinator.Coordinator()
    owner.add_event(simulator_event(sim1)).register_species(('A1', 'A2'))
    owner.add_event(ODEEvent(sim2, 0.01)).register_species(('E1', 'E2'))
    owner.initialize()

    data = []
    def log(owner):
        data.append((
            owner.t(),
            owner.get_value(Species("A1")),
            owner.get_value(Species("A2")),
            owner.get_value(Species("E1")),
            owner.get_value(Species("E2")),
            w2.get_value_exact(Species("A1")),
            ))

    log(owner)
    while owner.step(50):
        if owner.last_event.updated():
            log(owner)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt

    data = numpy.asarray(data)
    for i in range(1, 6):
        plt.plot(data.T[0], data.T[i], '-')
    # plt.show()
    plt.savefig('result.eps')
    numpy.savetxt("result.txt", data)


if __name__ == "__main__":
    # test1()
    test2()
