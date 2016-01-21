import libsbml


class SBMLDataSource(object):

    def __init__(self, filename):
        self.__data = libsbml.SBMLReader().readSBML(filename)
        self.__model = self.__data.getModel()

        # print(dir(self.__model))
        # print(self.__model.getListOfParameters())
        # print(self.__model.getListOfReactions())

    def initial_amounts(self):
        compartments = dict(self.compartments())

        for sp in self.__model.species:
            if sp.isSetInitialAmount():
                yield (sp.id, sp.initial_amount)
            elif sp.isSetInitialConcentration():
                yield (sp.id, sp.initial_concentration)
                # yield (sp.id, sp.initial_concentration * compartments[sp.compartment])

    def compartments(self):
        for comp in self.__model.compartments:
            yield (comp.id, comp.volume)

    def parameters(self):
        for p in self.__model.parameters:
            yield (p.id, p.value)

    def assignment_rules(self):
        for rule in self.__model.rules:
            if rule.isAssignment():
                yield (rule.variable, rule.formula)

    def reactions(self):
        for r in self.__model.reactions:
            reactants = [(reactant.species, reactant.stoichiometry)
                         for reactant in r.reactants]
            products = [(product.species, product.stoichiometry)
                        for product in r.products]

            formula = r.getKineticLaw().formula
            parameters = dict((p.id, p.value) for p in r.getKineticLaw().parameters)

            yield (reactants, products, formula, parameters)


if __name__ == '__main__':
    import sys
    # from operator import add
    # from functools import reduce, partial


    from ecell4 import *
    sbml = SBMLDataSource

    filename = sys.argv[1]

    y0 = dict(sbml(filename).initial_amounts())
    # y0.update(dict(sbml(filename).compartments()))
    print(y0)

    params = dict(sbml(filename).parameters())
    params.update(dict(sbml(filename).compartments()))
    params['compartment'] = 1.0
    print(params)

    with reaction_rules():
        # params.update(dict((var, _eval(formula)) for var, formula in sbml(filename).assignment_rules()))
        print(dict((var, _eval(formula)) for var, formula in sbml(filename).assignment_rules()))

        for reactants, products, formula, parameters in sbml(filename).reactions():
            parameters.update(params)

            (sum((_eval(sp) * coef for sp, coef in reactants), ~EmptySet)
                    > sum((_eval(sp) * coef for sp, coef in products), ~EmptySet) | _eval(formula, parameters))

            # _sum(reactants) > _sum(products) | _eval(formula, params)
            # _sum(_mul(reactants, reactant_coefficients)) > _sum(_mul(products, product_coefficients)) | _eval(formula, params)

    m = get_model()
    print([rr.as_string() for rr in m.reaction_rules()])

    run_simulation(60, model=m, y0=y0, opt_kwargs={'legend': False})
    # w = run_simulation(0.0032, model=m, y0=y0, species_list=['x1'], return_type='world')
    # y0 = dict((sp.serial(), w.get_value(sp)) for sp in w.list_species())
    # y0['k5'] = 1.55
    # run_simulation(0.1, model=m, y0=y0, species_list=['x1'])
