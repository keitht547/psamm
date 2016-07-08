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
import math
from itertools import count
from psamm.expression import boolean

from ..command import SolverCommandMixin, MetabolicMixin, Command, CommandError
from .. import fluxanalysis
from ..util import MaybeRelative
import csv
import math
import random
import psamm.lpsolver



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
        parser.add_argument('--transc-file', help='Enter path to transcriptomic data file',
        metavar='FILE')
        super(MadeFluxBalance, cls).init_parser(parser)


    def run(self):
        """Run MADE implementation."""
        x = self.parse_dict()
        var_dict1 = {}
        var_dict2 = {}
        gv1 = {}
        gv2 = {}
        trdict = {}
        linear_ineq_list = []
        var_gen = ('y{}'.format(i) for i in count(1))
        problem = self.flux_setup()
        for key, value in x.iteritems():
            e = boolean.Expression(value)
            exp_gene_string(e.base_tree(), var_gen, problem, var_dict1, var_dict2, key, gv1, gv2, trdict, linear_ineq_list)
            print (key,value) #Prints reaction ID and GPR associations
            print ' '
        print gv1
        print ' '
        print gv2

<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b
        thresh_v = self.minimum_flux()
        problem.prob.add_linear_constraints(thresh_v[1] >= thresh_v[0])
=======
        self.minimum_flux(problem)
>>>>>>> Improved commenting on numerous functions.

        if self._args.transc_file != None:
            gd = IDC(open_file(self))

        nat_model = self._model
        mm = nat_model.create_metabolic_model()
<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b

        info_dict = rxn_info(mm, problem)
        add_final_constraints(info_dict, problem, gv1, gv2)
        make_obj_fun(gv1, gv2, gd[2], gd[3], trdict, problem)
=======
        info_dict= rxn_info(mm, problem)
        make_obj_fun(reaction_dict, gv2, gd[2], gd[3], trdict)
        add_final_constraints(info_dict, problem, reaction_dict)
>>>>>>> Improved commenting on numerous functions.


    def parse_dict(self):
        ''' Using the model file,returns a dictionary with reactions as keys and
        their associated gene logic (i.e. (gene 1 and gene 2) or gene 3) as
        values of type str.'''
        gene_dict = {}
        for i in self._model.parse_reactions():
            gene_dict[i.id] = i.genes
        return gene_dict


    def flux_setup(self):
        '''Creates a flux balance problem using the model files in the current
        directory.'''
        model_path = self._model
        nat_model = self._model
        mm_model = nat_model.create_metabolic_model()
        solver = self._get_solver(integer=True)
        mm = mm_model.copy()
        p = fluxanalysis.FluxBalanceProblem(mm, solver)
        return p


    def minimum_flux(self, problem):
        '''Returns a biomass flux threshold that is a fraction of the maximum flux.
        Not in final working condition yet - notice te absence of a return.'''
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
<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b
        # print 'Ojective Reaction: {}'.format(obj_func)
        # print 'Objective flux: {}'.format(obj_flux)
        # print 'Biomass flux: {}'.format(Biomass)
        thresh_val = thresh.value*obj_flux
        return(thresh_val, obj_var)

def make_obj_fun(gv1, gv2, gp, gd, tr, problem):
    #gv1 = xi; gv2 = xi+1;
    #gp = gene probability (pvalue), gd = dictionary with increasing/decreasing expression values
=======
        print 'Ojective Reaction: {}'.format(obj_func)
        print 'Objective flux: {}'.format(obj_flux)
        print 'Biomass flux: {}'.format(Biomass)


def make_obj_fun(gv1, gv2, gp, gd, tr):
    '''Constructs the MADE objective funtion from dictionaries of LP variables.
    gv1 = xi, reaction_dict; gv2 = xi+1, gv2;
    gp = gene probability (pvalue), gd = dictionary with increasing/decreasing expression values'''
>>>>>>> Improved commenting on numerous functions.
    MILP_obj = 0.0
    for gene, var in tr.iteritems():
        print gene, var
        wp = gp[gene]
        print type(gp[gene])
        if gd[gene] == 1:
            MILP_obj = MILP_obj + (-math.log10(gp[gene]))*(gv2[var[1]] - gv1[var[0]])
        elif gd[gene] == -1:
            MILP_obj = MILP_obj + (-math.log10(gp[gene]))*(gv1[var[0]] - gv2[var[1]])
        elif gd[gene] == 0:
            if gv2[var[1]]- gv1[var[0]] <= 0:
                MILP_obj = MILP_obj + (-math.log10(gp[gene]))*(-(gv2[var[1]]-gv1[var[0]]))
            elif gv2[var[1]]- gv1[var[0]] >= 0:
                MILP_obj = MILP_obj + (-math.log10(gp[gene]))*(gv2[var[1]]-gv1[var[0]])
    print 'Objective Function: {}'.format(MILP_obj)
<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b
    problem.prob.define(('v', 'obj'))
    obj_var = problem.get_flux_var('obj')
    print obj_var
    problem.prob.add_linear_constraints(obj_var == MILP_obj)
    print obj_var
    problem.prob.add_linear_constraints(obj_var <= 1000)
    problem.maximize('obj')
    Biomass = problem.get_flux('obj')
    print 'Maxed Flux: {}'.format(Biomass)
    y = range(1, 12)
    for i in y:
        print(i)
        x_val = problem.get_ineq('y{}'.format(i))
        x2_val = problem.get_ineq('y{}.2'.format(i))
        print(x_val, x2_val)
    x = ['rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5']
    for j in x:
        print(j)
        print(problem.get_flux(j))
    #obj_var = p.get_flux_var('obj')

def exp_gene_string(A, var_gen, problem, var_dict1, var_dict2, name, gv1, gv2, tr, linear_ineq_list):
    '''Opens all containers, defines content, outputs the linear ineqs'''
    var_dict1[A] = name
=======

def exp_gene_string(A, var_gen, var_dict, name, linear_ineq_list, problem, reaction_dict, var_dict2, gv2, dict2):
    '''Opens all gene-logic containers, defines content, outputs the linear ineqs
    by calling bool_ineqs().  Sorts data into dictionaries that are used in other
    functions.  Is recursive. No output.'''
    var_dict[A] = name
>>>>>>> Improved commenting on numerous functions.
    problem.prob.define(("i", name))
    namevar = problem.get_ineq_var(name)
    gv1[name] = namevar

    var_dict2[A] = name +'.2'
    problem.prob.define(("i", name +'.2'))
    namevar2 = problem.get_ineq_var(name + '.2')
    gv2[name + '.2'] = namevar2

    if type(A) is boolean.Variable:
            print(A)
            str(A)
            tr.setdefault(A.symbol, []).append(name)
            tr.setdefault(A.symbol, []).append(name + '.2')
            print tr

    if type(A) is not boolean.Variable:
        exp_obj_name =var_dict1.get(A)
        children = []
        variable_names = []
        variable_names_gv2 = []
        for N,i in enumerate(A):
            children.append(i)
            newvar = next(var_gen)
            variable_names.append(newvar)
            variable_names_gv2.append(newvar + '.2')
            exp_gene_string(i, var_gen, problem, var_dict1, var_dict2, newvar, gv1, gv2, tr, linear_ineq_list)
            indent = (N+1) * '\t'
        for name1 in variable_names:
            problem.prob.define(("i",name1))
        for name2 in variable_names_gv2:
            problem.prob.define(("i", name2))

        print '{}Container Expression: '.format(indent), A
        print '{}Arguments: '.format(indent), children
        print '{}Variable Names: '.format(indent), variable_names

        if exp_obj_name is None:
            exp_obj_name = name
<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b
        x = bool_ineqs(A.cont_type(), A.contain(), variable_names, var_dict1, exp_obj_name, problem)
=======
        x = bool_ineqs(A.cont_type(), A.contain(), variable_names, var_dict,
        exp_obj_name, problem)
>>>>>>> Improved commenting on numerous functions.
        linear_ineq_list.append(x)

        print '{}Container Expression: '.format(indent), A
        print '{}Arguments: '.format(indent), children
        print '{}Variable Names: '.format(indent), variable_names_gv2
        y = bool_ineqs(A.cont_type(), A.contain(), variable_names_gv2, var_dict2,
        exp_obj_name, problem)
        linear_ineq_list.append(y)


def bool_ineqs(ctype, containing, names, dict_var, obj_name, problem):
    '''Input homogenous boolean.Expression (all ANDs or all ORs).  Adds the
    corresponding linear inequalities to the LP problem. No output.'''

    N = len(containing) # Length of the children list
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

    # Used for outputting a string of inequalities
    # ineq = [] #The list of inequalities to be returned
    # ineq1 = ' + '.join(x) # The first inequality
    # ineq.append(Y+relation1+ineq1+modify)
    # for j in range(N):
    #     if obj_name is not None:
    #         ineq.append(obj_name+relation2+x[j])
    #     else:
    #          ineq.append(obj_name+relation2+x[j])# Subsequent inequalities
    # print ineq

    yvar = problem.get_ineq_var(Y)
    RHSlist = []
    for name in names:
        xvar = problem.get_ineq_var(name)
        RHSlist.append(xvar)

    if label == 'or':
        or_group = None
        for ineq_var in RHSlist:
            linear_fxn(problem, 0 <= yvar <= 1)
            linear_fxn(problem, 0 <= ineq_var <= 1)
            #Individual variable inequalities
            orindiv = yvar >= ineq_var
            linear_fxn(problem, orindiv)
            print 'OR Constraint: ', orindiv
            #Container inequalities
            if or_group is None:
                or_group = ineq_var
            else:
                or_group = ineq_var + or_group
        or_cont = yvar <= or_group
        linear_fxn(problem, or_cont)
        print 'OR Constraint: ', or_cont

    if label == 'and':
        and_group = None
        for ineq_var in RHSlist:
            linear_fxn(problem, 0 <= yvar <= 1)
            linear_fxn(problem, 0 <= ineq_var <= 1)
            #Individual variable inequalities
            andindiv = yvar <= ineq_var
            linear_fxn(problem, andindiv)
            print 'AND Constraint: ', andindiv
            #Container inequalities
            if and_group is None:
                and_group = ineq_var
            else:
                and_group = ineq_var + and_group
        and_cont = yvar >= and_group - (N-1)
        linear_fxn(problem, and_cont)
        print 'AND Constraint: ', and_cont


def linear_fxn(lpp, linear_con):
    '''For adding linear inequalities as constraints in the LP problem specified.'''
    lpp.prob.add_linear_constraints(linear_con)

    # lpp = linear programming problem, linear_con = linear constraint


def flatten_list(biglist):
    '''Takes a list of lists and combines then into a singular list'''
    results = []
    for equations in biglist:
        for values in equations:
            results.append(values)
    return results

def open_file(self):
    '''Returns the contents of toy model file in a tuple of dictionaries.
    File Form: tsv format, FOUR Columns: (1) Gene name, (2) Condition 1 Data,
    (3) Condition 2 Data, (4) P-value of the transition 1->2.'''
    path = self._args.transc_file
    file1 = open(path)
    con1_dict = {}
    con2_dict = {}
    pval_dict = {}

    for row in csv.reader(file1, delimiter=str('\t')):
        try:
            con1_dict[row[0]] = float(row[1])
            con2_dict[row[0]] = float(row[2])
            pval_dict[row[0]] = float(row[3])
        except ValueError:
            print 'Cannot convert string to float',row[1], row[2], row[3]

    return con1_dict, con2_dict, pval_dict


def IDC(dicts, significance=0.05):
    '''Generates a dictionary with keys = genes and values = [-1, 0, +1]
    corresponding to significantly decreasing, constant, and inreasing expresson.
    P-values less than or equal to the significance are considered significant.'''
    con1 = dicts[0]
    con2 = dicts[1]
    pval = dicts[2]
    diff = {}
    for key in con1:
        if con2[key]-con1[key] == 0 or pval[key] > significance:
            diff[key] = 0
        else:
            diff[key] = int((con2[key]-con1[key])/abs(con2[key]-con1[key]))

    return con1,con2,pval,diff


def rxn_info(mm, problem):
    '''Takes a metabolic model and an LP Problem.
    Returns Dict:{rxn id: [low bound, high bound, fluxvar lp.Expression]}'''
    info = {}
    for rxn in mm.reactions:
<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b
        if mm.is_exchange(rxn) is False:
            info_list = []
            info_list.append(mm.limits.__getitem__(rxn).bounds[0])
            info_list.append(mm.limits.__getitem__(rxn).bounds[1])
            #   __getitem__ could be replaced by _create_bounds, or another function
            #   could be implemented.

            info_list.append(problem.get_flux_var(rxn))
            info[rxn] = info_list
    return info


def add_final_constraints(info_dict, problem, gv1, gv2):
=======
        info_list = []
        info_list.append(mm.limits.__getitem__(rxn).bounds[0])
        info_list.append(mm.limits.__getitem__(rxn).bounds[1])
        #   __getitem__ could be replaced by _create_bounds, or another function
        #   could be implemented.
        info_list.append(problem.get_flux_var(rxn))
        info[rxn] = info_list
    return info

def add_final_constraints(info_dict, problem, reaction_dict):
    '''Takes the output of rxn_info(), the LP problem, and the reaction dictionary.
    Adds contraints connecting flux variables, reactions, and flux bounds.'''
>>>>>>> Improved commenting on numerous functions.
    for rxn, info in info_dict.iteritems():
        print ''
        print rxn
        print info
<<<<<<< 8ff635122d8f6d48d50f7b0392fef03ebeb4c31b

        vmin = info[0]
        vmax = info[1]
        fluxvar = info[2]
        Y = gv1[rxn]
        # Z = gv2[rxn+'.2']
        print 'fluxvar: {}'.format(fluxvar)
        print 'min: {}'.format(vmin)
        print 'max: {}'.format(vmax)
        # print 'Z: {}'.format(Z)
        # except:
        #     print rxn, 'is not in the reaction dictionary.'
        linear_fxn(problem, fluxvar + Y*vmax >= 0)
        linear_fxn(problem, 0 <= fluxvar - Y*vmin)
        # linear_fxn(problem, fluxvar + Z*vmax >= 0)
        # linear_fxn(problem, 0 <= fluxvar - Z*vmin)
=======
        try:
            vmin = info[0]
            vmax = info[1]
            fluxvar = info[2]
            Y = reaction_dict[rxn]
            linear_fxn(problem, fluxvar + Y*vmax <= 0)
            linear_fxn(problem, 0 <= fluxvar - Y*vmin)
        except:
            print rxn, 'is not in the reaction dictionary.'
>>>>>>> Improved commenting on numerous functions.
