# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2016  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2016  Matthew Gentry <mgentry@umass.edu>

from __future__ import unicode_literals

import argparse
import time
import logging
import math
from itertools import count
from psamm.expression import boolean

from six import iteritems

from ..command import SolverCommandMixin, MetabolicMixin, Command, CommandError
from .. import fluxanalysis
from ..util import MaybeRelative
import csv
import psamm.lpsolver
from psamm.lpsolver import lp
from made import open_file




logger = logging.getLogger(__name__)


class eFluxBalance(MetabolicMixin, SolverCommandMixin, Command):
    """Run eflux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--loop-removal', help='Select type of loop removal constraints',
            choices=['none', 'tfba', 'l1min'], default='none')
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        parser.add_argument(
            '--flux-threshold',
            help='Enter minimum objective flux as a decimal or percent',
            type=MaybeRelative, default = MaybeRelative('90%'), nargs='?')
        parser.add_argument('--transc-file', help='Enter path to transcriptomic data file',
        metavar='FILE')
        super(eFluxBalance, cls).init_parser(parser)


    def run(self):
        """Run MADE implementation."""
        print ''
        print 'RUN STARTS HERE'
        print ''

        data = open_file(self)
        rxn_genelogic = self.gene_logic()
        minimum = self.minimum_flux()
        problem, mm =  self.flux_setup()
        ex_bounds = exchange_bounds(mm, problem)
        bounds = self.reaction_bounds(mm, problem)
        biomassflux = self.min_biomass(problem)

        rxn_exp, x = reaction_expression(rxn_genelogic, data)
        gene_bounds(mm, rxn_exp, x)
        bounds2 = self.reaction_bounds(mm, problem)
        biomass = self._model.get_biomass_reaction()
        solver = self._get_solver()
        problem2 = fluxanalysis.FluxBalanceProblem(mm, solver)
        problem2.maximize(biomass)

        print ''
        print 'Transcriptomic Data: ', data
        print ''
        print 'Reaction Gene Logic: ', rxn_genelogic
        print ''
        print 'Biomass threshold 90%: ', minimum
        print ''
        print 'Bounds:', bounds
        print ''
        #print 'Biomass: ', biomassflux
        print ''
        print 'A0, B0: ', ex_bounds
        print ''
        print 'Reaction Expression: ', rxn_exp
        print ''
        print 'New Bounds: '
        for rxn in mm.reactions:
            print rxn, mm.limits[rxn].bounds
        print ''
        print 'FBA:'
        for rxn in self._model.parse_reactions():
            print '{}\t{}\t{}\t{}'.format(rxn.id, problem2.get_flux(rxn.id),
            rxn.equation, rxn.genes)
        print 'Bounds 2: ', bounds2


    def gene_logic(self):
        ''' Using the model file,returns a dictionary with reactions as keys and
        their associated gene logic (i.e. (gene 1 and gene 2) or gene 3) as
        values of type str.'''
        gene_dict = {}
        for i in self._model.parse_reactions():
            gene_dict[i.id] = i.genes
        return gene_dict


    def reaction_bounds(self, mm, problem):
        '''Obtains reaction bounds, defines LP variables for all metabolic and
        exchange reactions, applies constraints, and returns a dictionary with
        all the bounds.'''
        bounds = {}
        for rxn in mm.reactions:
            bound = mm.limits[rxn].bounds
            problem.prob.define(('v', rxn))
            fluxvar = problem.get_flux_var(rxn)
            bounds[rxn] = bound
        return bounds



    def flux_setup(self):
        '''Creates a flux balance problem using the model files in the current
        directory.'''
        model_path = self._model
        nat_model = self._model
        mm_model = nat_model.create_metabolic_model()
        solver = self._get_solver(integer=True)
        mm = mm_model.copy()
        p = fluxanalysis.FluxBalanceProblem(mm, solver)
        return p, mm_model


    def minimum_flux(self):
        '''Returns a biomass flux threshold that is a fraction of the maximum flux.'''
        thresh = self._args.flux_threshold
        solver = self._get_solver(integer=True)
        mm_model = self._mm
        nat_model = self._model
        obj_func = nat_model.get_biomass_reaction()
        p = fluxanalysis.FluxBalanceProblem(mm_model, solver)
        p.maximize(obj_func)
        obj_flux = p.get_flux(obj_func)
        thresh_val = thresh.value*obj_flux
        return thresh_val


    def min_biomass(self, problem):
        '''Applies biomass minimum of 90% (defined in the default in init_parser)
        to the LP problem. Maximizes the biomass.'''
        biomass = self._model.get_biomass_reaction()
        threshold = self.minimum_flux()
        obj = problem.get_flux_var(biomass)
        problem.prob.add_linear_constraints(threshold <= obj)
        return biomass, threshold

def unpack(container, data):
    '''Recursively unpacks gene logic and returns the 'expression' of the
    reaction taken as an argument.  Takes a reaction boolean expression and
    the trancriptomic data.'''

    if isinstance(container, boolean.Variable):
        print str(container), data[str(container)]
        return data[str(container)]
    else:
        if isinstance(container, boolean.Or):
            x = []
            for i in container:
                x.append(unpack(i, data))
            print 'Or', sum(x)
            return sum(x)
        elif isinstance(container, boolean.And):
            x = []
            for i in container:
                x.append(unpack(i, data))
            print 'And', min(x)
            return min(x)


def exchange_bounds(mm, problem):
    '''Returns the baseline lower and upper bounds on the exchange reactions in
    the model.  Takes an instance of the metabolic model and the LP problem.'''
    A0 = {}
    B0 = {}
    w = {}
    exchange_rxns = []
    metabolic_rxns = []
    for rxn in mm.reactions:
        w[rxn] = 0
        if mm.is_exchange(rxn):
            exchange_rxns.append(rxn)
        else:
            metabolic_rxns.append(rxn)
    print 'exchange reactions: ',exchange_rxns
    print 'metabolic reactions: ',metabolic_rxns
    for exch in exchange_rxns:
        vmink = []
        vmaxk = []
        for rxn in metabolic_rxns:
            vmaxk.append(problem.flux_bound(rxn, 1))
            vmink.append(problem.flux_bound(rxn, -1))
        print exch, 'vmink -> vmaxk: ', vmink, vmaxk
        A0[exch] = min(vmink)
        B0[exch] = max(vmaxk)
    return A0, B0


def reaction_expression(gene_logic, data):
    '''Returns a reaction 'expression' dictionary of the form:
    {'con1 : {rxn1 : Xjl'}'}, where Xjl is the expression of reaction 1 under
    condition 1.

    gene-logic = dictionary connecting reactions to their gene logic
    data = the trancriptomic data from open_file().'''

    conditions = len(data)-1 #take one away to ignore the p-val column
    rxn_exp = {}
    x = []
    for con in range(conditions):
        key = str('con'+str(con+1))
        rxn_exp[key] = {}

        for rxn, logic in gene_logic.iteritems():
            exp = boolean.Expression(logic)
            Xjl = unpack(exp.base_tree(), data[con])
            rxn_exp[key][rxn] = Xjl
            x.append(Xjl)
            print 'Xjl =', Xjl
            print ''
    return rxn_exp, x


def gene_bounds(mm, rxn_exp, x):
    '''Alters the bounds on metabolic reactions by a factor equal to the reaction
    expression, relative to the maximum expression accross all conditions.
    Takes the metabolic model and the output of reaction_expression.  x is a
    list of all expressions over all reactions and all conditions.'''
    maxx = max(x)
    rxns = rxn_exp['con1']
    for rxn, Xjl in rxns.iteritems():
        print rxn, Xjl
        factor = Xjl/maxx
        a0 = mm.limits[rxn].lower
        mm.limits[rxn].lower = a0*factor
        b0 = mm.limits[rxn].upper
        mm.limits[rxn].upper = b0*factor
