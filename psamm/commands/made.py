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

import argparse
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
        linear_ineq_list = []
        var_gen = ('y{}'.format(i) for i in count(1))
        problem = self.flux_setup()
        for key, value in x.iteritems():
            e = boolean.Expression(value)
            exp_gene_string(e.base_tree(), var_gen, var_dict, key, linear_ineq_list, problem)
            print (key,value) #Prints reaction ID and gene string
            print ' '
        self.minimum_flux()

        # master_ineq_list = flatten_list(linear_ineq_list) #Complete list of inequalities
        # print master_ineq_list



    def parse_dict(self):
        gene_dict = {}
        for i in self._model.parse_reactions():
            gene_dict[i.id] = i.genes
        return(gene_dict)


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
        
        obj_var = p.get_flux_var(obj_func)
        linear_fxn(p, obj_var >= thresh.value*obj_flux)
        p.minimize_l1()
        Biomass = p.get_flux(obj_func)

        print obj_func # reaction
        print obj_flux #after maximized
        print thresh.value*obj_flux #maxed flux * threshold
        print Biomass # after minimzed
        # return obj_func
        # return obj_flux
        # return thresh.value*obj_flux


    def flux_setup(self):
        '''Creates a flux balance problem'''
        model_path = self._model
        nat_model = self._model
        mm_model = nat_model.create_metabolic_model()
        solver = self._get_solver(integer=True)
        mm = mm_model.copy()
        p = fluxanalysis.FluxBalanceProblem(mm, solver)
        return p


def exp_gene_string(A, var_gen, var_dict, name, linear_ineq_list, problem):
    '''Opens and identifies all containers with the variables and arguments
        as well outputs the associated inequalities'''
    var_dict[A] = name
    problem.prob.define(("i", name))
    if type(A) is not boolean.Variable:
        exp_obj_name =var_dict.get(A)
        children = []
        variable_names = []
        for N,i in enumerate(A):
            children.append(i)
            q = next(var_gen)
            variable_names.append(q)
            exp_gene_string(i, var_gen, var_dict, q, linear_ineq_list, problem)
            indent = (N+1) * '\t'
        for j in variable_names:
            problem.prob.define(("i",j))
        # for i in problem.prob._variables:
        #     print(i)

        if i in variable_names:
             print '{}Var Name: '.format(indent),variable_names(i)
        print '{}Container Expression: '.format(indent), A
        print '{}Arguments: '.format(indent), children
        print '{}Variable Names: '.format(indent), variable_names

        if exp_obj_name is None:
            exp_obj_name = name
        x = bool_ineqs(A.cont_type(), A.contain(), variable_names, var_dict, exp_obj_name, problem) #Prints the inequalities in list form
        linear_ineq_list.append(x)



def bool_ineqs(ctype, containing, names, dict_var, obj_name, problem):
    '''Input homogenous boolean.Expression type.
    Returns a list of corresponding unicode inequalities'''

    N = len(containing) # Length of the chilren list
    if isinstance(ctype, boolean.And):
        label = 'and'
        relation1 = ' >= '
        relation2 = ' <= '
        modify = ' - '+unicode(N-1) # one less than the number of ands or ors

    elif isinstance(ctype, boolean.Or):
        label = 'or'
        relation1 = ' <= '
        relation2 = ' >= '
        modify = ''
    elif isinstance(ctype, boolean.Variable):
        raise ValueError('Argument contains only variables, no operators')

    x = names # A list of the unicode characters of the variables in the expression
    if ctype in dict_var.keys():
        Y = dict_var[ctype]
    else:
        Y = 'Y'

    ineq = [] #The list of inequalities to be returned
    ineq1 = ' + '.join(x) # The first inequality
    ineq.append(Y+relation1+ineq1+modify)
    for j in range(N):
        if obj_name is not None:
            ineq.append(obj_name+relation2+x[j])
        else:
             ineq.append(obj_name+relation2+x[j])# Subsequent inequalities
    # return ineq


    '''The following was appended on 6/24/16'''
    yvar = problem.get_ineq_var(Y)
    rightsidelist = []
    for i in names:
        xvar = problem.get_ineq_var(i)
        rightsidelist.append(xvar)
    # print yvar
    # print rightsidelist

    if label == 'or':
        M = None
        for k in rightsidelist:
            L = yvar >= k
            linear_fxn(problem, L)
            print L
            if M is None:
                M = k
            else:
                M = k+M
        Q = yvar <= M
        linear_fxn(problem, Q)
        print Q

    if label == 'and':
        S = None
        for k in rightsidelist:
            R = yvar <= k
            linear_fxn(problem, R)
            print R
            if S is None:
                S = k
            else:
                S = k + S
        T = yvar >= S - (N-1)
        linear_fxn(problem, T)
        print T


def linear_fxn(lpp, linear_con):
    '''Created for ease of adding linear constraints'''
    # lpp = linear programming or 'p' in this case, linear_con = linear constraint
    lpp.prob.add_linear_constraints(linear_con)


def flatten_list(biglist):
    '''Takes a list of lists and combines then into a singular list'''
    results = []
    for equations in biglist:
        for values in equations:
            results.append(values)
    return results
