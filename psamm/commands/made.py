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
        parser = self.parse_dict()
        gene_var1 = {} #Dictionary; key:value = gene expression objects:their new variable id, first set
        gene_var2 = {} #Dictionary; key:value = gene expression objects:their new variable id, second set
        var_ineqvar1 = {} #Dictionary; key:value = new variable ids:their defined inequality variable, first set
        var_ineqvar2 = {} #Dictionary; key:value = new variable ids:their defined inequality variable, second set
        gvdict = {} #Dictionary; key:value = gene_id:defined variable ids from both sets (each key has 2 values)

        var_gen = ('y{}'.format(i) for i in count(1))
        problem = self.flux_setup()
        for rxn_id, GPR in parser.iteritems():
            e = boolean.Expression(GPR)
            exp_gene_string(e.base_tree(), var_gen, problem, rxn_id, gene_var1, gene_var2, var_ineqvar1, var_ineqvar2, gvdict)
            print (rxn_id,GPR)
            print ' '
        print var_ineqvar1
        print ' '
        print var_ineqvar2


        thresh_v = self.minimum_flux()
        problem.prob.add_linear_constraints(thresh_v[1] >= thresh_v[0])

        if self._args.transc_file != None:
            gene_data = IDC(open_file(self))

        nat_model = self._model
        mm = nat_model.create_metabolic_model()


        info_dict = rxn_info(mm, problem)
        add_final_constraints(info_dict, problem, var_ineqvar1, var_ineqvar2)
        make_obj_fun(var_ineqvar1, var_ineqvar2, gene_data[2], gene_data[3], gvdict, problem)





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


    def minimum_flux(self):
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

        # print 'Ojective Reaction: {}'.format(obj_func)
        # print 'Objective flux: {}'.format(obj_flux)
        # print 'Biomass flux: {}'.format(Biomass)
        thresh_val = thresh.value*obj_flux
        return(thresh_val, obj_var)

def make_obj_fun(var_ineqvar1, var_ineqvar2, gene_pval, gene_data, gvdict, problem):
    '''Constructs the MADE objective funtion from dictionaries of LP variables.
    var_ineqvar1 = xi, ; var_ineqvar2 = xi+1, var_ineqvar2; gene_pval = gene probability (pvalue),
    gene_data = dictionary with increasing/decreasing expression values'''
    MILP_obj = 0.0
    # test = 0.0
    for gene, var in gvdict.iteritems():
        print gene, var
        wp = gene_pval[gene]
        print type(gene_pval[gene])
        if gene_data[gene] == 1:
            MILP_obj = MILP_obj + (-math.log10(gene_pval[gene]))*(var_ineqvar2[var[1]] - var_ineqvar1[var[0]])
        elif gene_data[gene] == -1:
            MILP_obj = MILP_obj + (-math.log10(gene_pval[gene]))*(var_ineqvar1[var[0]] - var_ineqvar2[var[1]])
        elif gene_data[gene] == 0:
            if var_ineqvar2[var[1]]- var_ineqvar1[var[0]] <= 0:
                MILP_obj = MILP_obj - (-math.log10(gene_pval[gene]))*(-(var_ineqvar2[var[1]]-var_ineqvar1[var[0]]))
            elif var_ineqvar2[var[1]]- var_ineqvar1[var[0]] >= 0:
                MILP_obj = MILP_obj - (-math.log10(gene_pval[gene]))*(var_ineqvar2[var[1]]-var_ineqvar1[var[0]])
    print 'Objective Function: {}'.format(MILP_obj)

    y = range(1, 12)
    for i in y:
        y_val = problem.get_ineq_var('y{}'.format(i))
        y2_val = problem.get_ineq_var('y{}.2'.format(i))

        linear_fxn(problem, y_val <= 1)
        linear_fxn(problem, 0 <= y_val)
        linear_fxn(problem, y2_val <= 1)
        linear_fxn(problem, 0 <= y2_val)


    problem.prob.define(('v', 'obj'))
    obj_var = problem.get_flux_var('obj')
    problem.prob.add_linear_constraints(obj_var == MILP_obj)
    # print obj_var
    problem.prob.add_linear_constraints(obj_var <= 50000)
    problem.maximize('obj')
    obj_flux = problem.get_flux('obj')

    print 'Maxed Flux: {}'.format(obj_flux)
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


def exp_gene_string(exp_obj, var_gen, problem, new_var_id, gene_var1, gene_var2, var_ineqvar1, var_ineqvar2, gvdict):
    '''Opens all gene-logic containers, defines content, outputs the linear ineqs
    by calling bool_ineqs().  Sorts data into dictionaries that are used in other
    functions.  Is recursive. No output.'''

    gene_var1[exp_obj] = new_var_id
    problem.prob.define(("i", new_var_id))
    new_var_ineq1 = problem.get_ineq_var(new_var_id)
    var_ineqvar1[new_var_id] = new_var_ineq1
    0<= var_ineqvar1[new_var_id] <=1
    linear_fxn(problem, 0<= var_ineqvar1[new_var_id] <=1)

    gene_var2[exp_obj] = new_var_id +'.2'
    problem.prob.define(("i", new_var_id +'.2'))
    new_var_ineq2 = problem.get_ineq_var(new_var_id + '.2')
    var_ineqvar2[new_var_id + '.2'] = new_var_ineq2
    0<= var_ineqvar2[new_var_id+ '.2'] <=1
    linear_fxn(problem, 0<= var_ineqvar2[new_var_id+ '.2'] <=1)

    if type(exp_obj) is boolean.Variable:
        print(exp_obj)
        str(exp_obj)
        gvdict.setdefault(exp_obj.symbol, []).append(new_var_id)
        gvdict.setdefault(exp_obj.symbol, []).append(new_var_id + '.2')
        print gvdict

    if type(exp_obj) is not boolean.Variable:
        exp_obj_name =gene_var1.get(exp_obj)
        arguments = []
        new_var_names1 = []
        new_var_names2 = []
        for N,i in enumerate(exp_obj):
            arguments.append(i)
            newvar = next(var_gen)
            new_var_names1.append(newvar)
            new_var_names2.append(newvar + '.2')
            exp_gene_string(i, var_gen, problem, newvar, gene_var1, gene_var2, var_ineqvar1, var_ineqvar2, gvdict)
            indent = (N+1) * '\t'
        for name1 in new_var_names1:
            problem.prob.define(("i",name1))
        for name2 in new_var_names2:
            problem.prob.define(("i", name2))

        print '{}Container Expression: '.format(indent), exp_obj
        print '{}Arguments: '.format(indent), arguments
        print '{}Variable Names: '.format(indent), new_var_names1

        if exp_obj_name is None:
            exp_obj_name = new_var_id

        x = bool_ineqs(exp_obj.cont_type(), exp_obj.contain(), new_var_names1, gene_var1, exp_obj_name, problem)


        print '{}Container Expression: '.format(indent), exp_obj
        print '{}Arguments: '.format(indent), arguments
        print '{}Variable Names: '.format(indent), new_var_names2
        y = bool_ineqs(exp_obj.cont_type(), exp_obj.contain(), new_var_names2, gene_var2,
        exp_obj_name, problem)


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
    print 'nimei: {}'.format(linear_con)
    lpp.prob.add_linear_constraints(linear_con)
    # lpp = linear programming problem, linear_con = linear constraint


def open_file(self):
    '''Returns the contents of model file in a tuple of dictionaries.
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

        if mm.is_exchange(rxn) is False:
            info_list = []
            info_list.append(mm.limits.__getitem__(rxn).bounds[0])
            info_list.append(mm.limits.__getitem__(rxn).bounds[1])
            #   __getitem__ could be replaced by _create_bounds, or another function
            #   could be implemented.
            info_list.append(problem.get_flux_var(rxn))
            info[rxn] = info_list
    return info


def add_final_constraints(info_dict, problem, var_ineqvar1, var_ineqvar2):
    '''Takes the output of rxn_info, the LP Problem, and the binary dictionaries
    of each condition.  Adds constraints connecting flux variables, reactions,
    and their flux bounds.'''
    for rxn, info in info_dict.iteritems():
        print ''
        vmin = info[0]
        vmax = info[1]
        fluxvar = info[2]
        Y = var_ineqvar1[rxn]
        Z = var_ineqvar2[rxn+'.2']
        print 'fluxvar: {}'.format(fluxvar)
        print 'min: {}'.format(vmin)
        print 'max: {}'.format(vmax)
        print 'Z: {}'.format(Z)

        linear_fxn(problem, fluxvar + (1-Y)*vmax <= vmax)
        linear_fxn(problem, fluxvar + (1-Y)*vmin >= vmin)
        linear_fxn(problem, fluxvar + (1-Z)*vmax <= vmax)
        linear_fxn(problem, fluxvar + (1-Z)*vmin >= vmin)
