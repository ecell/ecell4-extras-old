import libsbml


class SBMLDataSource(object):

    def __init__(self, filename):
        self.__data = libsbml.SBMLReader().readSBML(filename)
        self.__model = self.__data.getModel()

        # print(dir(self.__model))
        # print(self.__model.getListOfParameters())
        # print(self.__model.getListOfReactions())

    def species(self):
        for sp in self.__model.species:
            if sp.isSetInitialAmount():
                yield (sp.id, sp.initial_amount)
            # if sp.isSetInitialConcentration():
            #     yield (sp.id, sp.initial_amount, sp.initial_concentration)

    def parameters(self):
        for p in self.__model.parameters:
            yield (p.id, p.value)

    def reactions(self):
        for r in self.__model.reactions:
            reactants = [(reactant.species, reactant.stoichiometry)
                         for reactant in r.reactants]
            products = [(product.species, product.stoichiometry)
                        for product in r.products]
            formula = r.getKineticLaw().formula
            yield (reactants, products, formula)


if __name__ == '__main__':
    import sys
    # from operator import add
    # from functools import reduce, partial


    from ecell4 import *
    sbml = SBMLDataSource

    filename = sys.argv[1]

    y0 = dict(sbml(filename).species())
    print(y0)
    params = dict(sbml(filename).parameters())
    print(params)

    with reaction_rules():
        for reactants, products, formula in sbml(filename).reactions():
            # print(reactants, products, formula)
            # concat(reactants) > concat(products) | _eval(formula)

            (sum((_eval(sp) * coef for sp, coef in reactants), ~X)
                    > sum((_eval(sp) * coef for sp, coef in products), ~X) | _eval(formula, params))

            # line = ''
            # if len(reactants) > 0:
            #     line += '+'.join('{}*{}'.format(sp, coef) for sp, coef in reactants)
            # else:
            #     line += '~X'
            # if len(products) > 0:
            #     line += '>' + '+'.join('{}*{}'.format(sp, coef) for sp, coef in products)
            # else:
            #     line += '>~X'
            # line += '|' + formula
            # evaluate(line)

    m = get_model()
    print([rr.as_string() for rr in m.reaction_rules()])

    # with reaction_rules():
    #     for reactants, products, formula in sbml(filename).reactions():
    #         # print(reactants, products, formula)
    #         print(repr(evaluate(formula)))
