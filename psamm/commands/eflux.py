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
        bounds = self.reaction_bounds(mm, problem)
        biomassflux = self.min_biomass(problem)
        print ''
        for rxn, logic in rxn_genelogic.iteritems():
            exp = boolean.Expression(logic)
            print rxn
            Xlj = unpack(exp.base_tree(), data[0])
            print 'Xlj = ', Xlj
            print ''

        print ''
        print 'Transcriptomic Data: ', data
        print ''
        print 'Reaction Gene Logic: ', rxn_genelogic
        print ''
        print 'Biomass threshold 90%: ', minimum
        print ''
        print 'Bounds:', bounds
        print ''
        print 'Biomass: ', biomassflux
        print ''
        exchange_bounds(mm, problem)


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
        exchange reactions, and returns a dictionary with all the bounds.'''
        bounds = {}
        for rxn in mm.reactions:
            bounds[rxn] = rxn_bounds(mm, rxn)
            problem.prob.define(('v', rxn)) #defines LP variables for all reactions
            fluxvar = problem.get_flux_var(rxn)
            problem.prob.add_linear_constraints(fluxvar <= bounds[rxn][1], fluxvar >= bounds[rxn][0])
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
        to the LP problem. Applies arbitrary cap of 2000. Maximizes the biomass.'''
        biomass = self._model.get_biomass_reaction()
        threshold = self.minimum_flux()
        problem.prob.define(('v', biomass))
        obj = problem.get_flux_var(biomass)
        problem.prob.add_linear_constraints(obj <= 2000)
        problem.prob.add_linear_constraints(threshold <= obj)
        problem.maximize(biomass)
        return biomass, problem.get_flux(biomass)


def unpack(container, data):
    '''Recursively unpacks gene logic and returns the 'expression' of the
    reaction taken as an argument.  Takes a reaction boolean expression and
    the trancriptomic data.'''

    if isinstance(container, boolean.Variable):
        print str(container)
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


def rxn_bounds(mm, rxn):
    '''Takes a metabolic model, a reaction and an LP Problem.
    Returns (low bound, high bound)'''

    lower_bound = mm.limits._create_bounds(rxn).bounds[0]
    upper_bound = mm.limits._create_bounds(rxn).bounds[1]
        #   __getitem__ could be replaced by _create_bounds, or another function
        #   could be implemented.
    return lower_bound, upper_bound

def exchange_bounds(mm, problem):
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
            problem.maximize(rxn)
            vmaxk.append(problem.get_flux(rxn))
            w[rxn] = 1
            problem.minimize_l1(w)
            vmink.append(problem.get_flux(rxn))
            w[rxn] = 0
        print exch, 'vmink -> vmaxk: ', vmink, vmaxk
        A0[exch] = min(vmink)
        B0[exch] = max(vmaxk)
    print 'A0: ', A0
    print 'B0: ', B0
