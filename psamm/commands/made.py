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

from __future__ import unicode_literals

import time
import logging
from itertools import count
from psamm.expression import boolean

from ..command import SolverCommandMixin, MetabolicMixin, Command, CommandError
from .. import fluxanalysis
from ..util import MaybeRelative


logger = logging.getLogger(__name__)


class MadeFluxBalance(MetabolicMixin, SolverCommandMixin, Command):
    """Run MADE flux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--loop-removal', help='Select type of loop removal constraints',
            choices=['none', 'tfba', 'l1min'], default='none')
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        parser.add_argument(
            '--flux-threshold',
            help='Enter maximum objective flux as a decimal or percent',
            type=MaybeRelative, default = MaybeRelative('100%'), nargs='?')
        super(MadeFluxBalance, cls).init_parser(parser)

    def run(self):
        """Run MADE implementation."""
#Reaction string
        x = self.parse_dict()
        var_dict = {}
        var_gen = ('y{}'.format(i) for i in count(1))
        for key, value in x.iteritems():
            e = boolean.Expression(value)
            exp_gene_string(e.base_tree(), var_gen, var_dict, key)
            print (key,value) #Prints reaction ID and gene string
            print ' '
            print ' '


    def parse_dict(self):
        gene_dict = {}
        for i in self._model.parse_reactions():
            gene_dict[i.id] = i.genes
        return(gene_dict)


def minimum_flux(self):
    '''Returns a biomass flux threshold that is a fraction of the maximum flux.'''
    q = self._args.flux_threshold
    solver = self._get_solver(integer=True)
    mm_model = self._mm
    nat_model = self._model
    obj_func = nat_model.get_biomass_reaction()
    p = fluxanalysis.FluxBalanceProblem(mm_model, solver)
    p.maximize(obj_func)
    flux = p.get_flux(obj_func)
    return q.value*flux


def bool_ineqs(exp):
    '''Input homogenous boolean.Expression type.
    Returns a list of corresponding unicode inequalities'''
    #print exp
    N = len(exp[1]) # Length of the chilren list
    if isinstance(exp[0], boolean.And):
        label = 'and'
        relation1 = ' >= '
        relation2 = ' <= '
        modify = ' - '+unicode(N-1) # one less than the number of ands or ors

    elif isinstance(exp[0], boolean.Or):
        label = 'or'
        relation1 = ' <= '
        relation2 = ' >= '
        modify = ''
    elif isinstance(exp[0], boolean.Variable):
        raise ValueError('Argument contains only variables, no operators')

    x = exp[2] # A list of the unicode characters of the variables in the expression
    ineq = [] #The list of inequalities to be returned
    ineq1 = ' + '.join(x) # The first inequality
    if exp[0] in exp[3].keys():
        Y = exp[3][exp[0]]
    else:
        Y = 'Y'
    ineq.append(Y+relation1+ineq1+modify)
    for j in range(N):
        if exp[4] is not None:
            ineq.append(exp[4]+relation2+x[j])
        else:
             ineq.append(exp[4]+relation2+x[j])# Subsequent inequalities
    return ineq



def exp_gene_string(A, var_gen, var_dict, name):
    var_dict[A] = name

    if type(A) is not boolean.Variable:
        exp_obj_name =var_dict.get(A)
        children = []
        variable_names = []
        for N,i in enumerate(A):
            children.append(i)
            q = next(var_gen)
            variable_names.append(q)
            exp_gene_string(i, var_gen, var_dict, q)
            indent = (N+1) * '\t'
        if i in variable_names:
             print '{}Var Name: '.format(indent),variable_names(i)

        print '{}Container Expression: '.format(indent), A
        print '{}Arguments: '.format(indent), children
        print '{}Variable Names: '.format(indent), variable_names

        if exp_obj_name is None:
            exp_obj_name = name
        Expression = [A.cont_type(), A.contain(), variable_names, var_dict, exp_obj_name]
        print bool_ineqs(Expression) #Prints the inequalities. List form
