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
        exp_and = boolean.Expression('A and B and C and D')
        exp_or = boolean.Expression('E or F or G or H')
        print bool_ineqs(exp_or)

#Reaction string
        x = self.parse_dict()
        var_gen = ('y{}'.format(i) for i in count(1))
        for key, value in x.iteritems():
            e = boolean.Expression(value)
            exp_gene_string(e.base_tree(), var_gen)
            print (key,value)
        # for key, value in x.iteritems():
        #     get_root(value)

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


def get_root(G):
    e = boolean.Expression(G)
    #print(e.base_tree())
    if type(e.base_tree()) is boolean.And:
        and_list = []
        for i in e.base_tree().contain():
            and_list.append(i)
        y1 = None
        for j in and_list:
            print('y1 <= {}'.format(j))
            if y1 is None:
                y1 = str(j)+' - {}'.format(len(and_list)-1)
            else:
                y1 = '{} + '.format(str(j))+y1
        y1 = 'y1 >= {}'.format(y1)
        print(y1)


def bool_ineqs(exp):
    '''Input homogenous boolean.Expression type.
    Returns a list of corresponding unicode inequalities'''

    N = len(exp.variables) # Counts the number of varables
    if isinstance(exp.base_tree(), boolean.And):
        label = 'and'
        relation1 = ' >= '
        relation2 = ' <= '
        modify = ' - '+unicode(N-1) # one less than the number of ands or ors

    elif isinstance(exp.base_tree(), boolean.Or):
        label = 'or'
        relation1 = ' <= '
        relation2 = ' >= '
        modify = ''
    elif isinstance(exp.base_tree(), boolean.Variable):
        raise ValueError('Argument contains only variables, no operators')

    x = [] # A list of the unicode characters of the variables in the expression
    for i in exp.base_tree().contain():
        x.append(i.symbol)

    ineq = [] #The list of inequalities to be returned
    ineq1 = ' + '.join(x) # The first inequality
    ineq.append('Y'+relation1+ineq1+modify)
    for j in range(N):
        ineq.append('Y'+relation2+x[j]) # Subsequent inequalities

    return ineq

def exp_gene_string(A, var_gen):
    if type(A) is not boolean.Variable:
        print type(A)
        children = []
        variable_names = []
        for i in A:
            children.append(i)
            variable_names.append(next(var_gen))
            exp_gene_string(i, var_gen)
            if i in variable_names:
                 print variable_names(i)
            print(children)
            print(variable_names)

    else:
        print A


    #
    # roots = []
    # containers = []
    # roots.append(e._root)
    # containers.append(e._root.contain())

    # print roots
    # # for index, item in enumerate(roots):
    # #     print index+1, item
    # print containers

    #if term in e._root not boolean.Or:
    #print e._root
    #print e.base_tree()
    # if type(e.base_tree()) is not boolean.Variable:
    #     for i in e.base_tree():
    #         if type(i) is not boolean.Variable:
    #             for j in i.contain():
    #                 if type(j) is not boolean.Variable:
    #                     for k in j.contain():
    #                         print(k)
    #                 print(j)
    #         else:
    # #             print(i)
    #         if type(i) is boolean.Or:
    #             print "or"
    #         elif type(i) is boolean.And:
    #             print "and"
    #         else:
    #             print "variable"
    # else:
    #     print e.base_tree()
