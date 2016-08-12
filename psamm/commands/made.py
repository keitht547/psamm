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

'''Implementation of Metabolic Adjustment by Differential Expression (MADE)'''

from __future__ import unicode_literals
import sys
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
from psamm.lpsolver import lp

# Module-level logging
logger = logging.getLogger(__name__)


class MadeFluxBalance(MetabolicMixin, SolverCommandMixin, Command):
    """Run MADE flux balance analysis on the model.

    Args:
        gene_var1 = Dictionary, key:value = gene expression objects:their new variable id, first set
        gene_var2 = Dictionary; key:value = gene expression objects:their new variable id, second set
        var_ineqvar1 = xi; Dictionary, key:value = new variable ids:their defined inequality variable, first set
        var_ineqvar2 = xi+1; Dictionary, key:value = new variable ids:their defined inequality variable, second set
        gene_pval = Dictionary, key:value = gene ID:gene fold change probability (pvalue)
        gene_diff = Dictionary, key:value = gene ID: binary up/down/constant regulation values
        gvdict = Dictionary, key:value = gene ID:defined variable ids from both sets (each key has 2 values)
        problem = Flux balance problem
    """

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
        parser.add_argument(
            '--fva', help='Enable FVA',
            action='store_true')
        super(MadeFluxBalance, cls).init_parser(parser)



    def run(self):
        """Run MADE implementation."""
        parser = self.parse_dict()
        gene_var1 = {}
        gene_var2 = {}
        var_ineqvar1 = {}
        var_ineqvar2 = {}
        gvdict = {}

        nat_model = self._model
        biomass_fun = nat_model.get_biomass_reaction()
        var_gen = ('y{}'.format(i) for i in count(1))
        problem_mm = self.flux_setup()
        start_time = time.time()
        mm = problem_mm[1]
        problem = problem_mm[0]
        thresh_v = self.minimum_flux()
        linear_constraints(problem, problem.get_flux_var(biomass_fun) >= thresh_v)

        for rxn_id, GPR in parser.iteritems():
            e = boolean.Expression(GPR)
            exp_gene_string(e.base_tree(), var_gen, problem, rxn_id, gene_var1, gene_var2, var_ineqvar1, var_ineqvar2, gvdict)


        if self._args.transc_file != None:
            gene_data = IDC(open_file(self))
        nat_model = self._model

        info_dict = rxn_info(mm, problem)
        add_final_constraints(info_dict, problem, var_ineqvar1, var_ineqvar2)
        final = make_obj_fun(var_ineqvar1, var_ineqvar2, gene_data[2], gene_data[3], gvdict, problem, biomass_fun)

        biomass_function = nat_model.get_biomass_reaction()
        for reaction in nat_model.parse_reactions():
            if reaction.id in nat_model.parse_model():
                if self._args.fva is True:
                    flux = final.get_flux_var(biomass_function)
                    linear_constraints(final, flux >= 0.99*final.get_flux(biomass_function))
                    print('{}\t{}\t{}\t{}\t{}'.format(reaction.id, final.flux_bound(reaction.id, -1),
                    final.flux_bound(reaction.id, 1), str(reaction.equation), reaction.genes))
                else:
                    print('{}\t{}\t{}\t{}'.format(reaction.id, final.get_flux(reaction.id), str(reaction.equation), reaction.genes))

        print('Objective Flux: {}'.format(final.get_flux(biomass_function)))
        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))


    def parse_dict(self):
        ''' Using the reaction file called inside of the model file, it returns a dictionary with reaction IDs as keys and
        their associated gene-protein reaction (GPR) logic (i.e. (gene 1 and gene 2) or gene 3) as
        values of type str.'''
        gene_dict = {}
        for i in self._model.parse_reactions():
            if i.genes is not None:
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


def linear_constraints(LPproblem, linear_constraint):
    '''For easily adding linear constraints to the linear programming problem.'''
    LPproblem.prob.add_linear_constraints(linear_constraint)

def make_obj_fun(var_ineqvar1, var_ineqvar2, gene_pval, gene_diff, gvdict, problem, biomass_fun):
    '''Constructs the MADE objective funtion from dictionaries of LP variables.

    Objective function consists of the summation of three functions dependent on the
    up/down regulation of gene expression between conditions. The functions contain
    a weighting function, and the difference between the binary representations of
    condition 1 and condition 2.
    '''

    MILP_obj = 0.0
    I = 0.0 # Increasing gene expression
    D = 0.0 # Decreasing gene expression
    C = 0.0 # Constant gene expression

    for gene, var in gvdict.iteritems():
        if gene in gene_pval.keys():
            wp = gene_pval[gene]
            if wp <= 2.2204460492e-16:
                wp = 2.2204460492e-16 #Limitation of math.log()
            else:
                wp = wp
            if gene_diff[gene] == 1:
                I = I + (-math.log10(wp))*(var_ineqvar2[var[1]] - var_ineqvar1[var[0]])
            elif gene_diff[gene] == -1:
                D = D + (-math.log10(wp))*(var_ineqvar1[var[0]] - var_ineqvar2[var[1]])
            elif gene_diff[gene] == 0:
                C = C + (-math.log10(wp))*((var_ineqvar2[var[1]]-var_ineqvar1[var[0]]))

            MILP_obj = I + D - C

    problem.prob.define(('v', 'object'))
    obj_var = problem.get_flux_var('object')
    linear_constraints(problem, obj_var == MILP_obj)
    problem.maximize('object')
    obj_flux = problem.get_flux('object')
    print obj_flux
    linear_constraints(problem, obj_var == obj_flux)
    problem.maximize(biomass_fun) # Must be changed based on ID of desired reaction
    return(problem)


def exp_gene_string(exp_obj, var_gen, problem, new_var_id, gene_var1, gene_var2, var_ineqvar1, var_ineqvar2, gvdict):
    '''Opens all gene-logic containers, defines content, outputs the linear inequalities
    by calling bool_ineqs().  Sorts data into dictionaries that are used in other
    functions.  Is recursive. No output.

    Args:
        exp_obj: All of the expression objects (genes, AND, OR)
        var_gen: Counter used for relabeling the genes and arguments as variables
        new_var_id: Variable ID, also includes original reaction ID for first layer
    '''

    gene_var1[exp_obj] = new_var_id
    problem.prob.define(("i", new_var_id), types = lp.VariableType.Binary)
    new_var_ineq1 = problem.get_ineq_var(new_var_id)
    var_ineqvar1[new_var_id] = new_var_ineq1

    gene_var2[exp_obj] = new_var_id +'.2'
    problem.prob.define(("i", new_var_id +'.2'), types = lp.VariableType.Binary)
    new_var_ineq2 = problem.get_ineq_var(new_var_id + '.2')
    var_ineqvar2[new_var_id + '.2'] = new_var_ineq2


    if type(exp_obj) is boolean.Variable:
        str(exp_obj)
        gvdict.setdefault(exp_obj.symbol, []).append(new_var_id)
        gvdict.setdefault(exp_obj.symbol, []).append(new_var_id + '.2')


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
            problem.prob.define(("i",name1),types = lp.VariableType.Binary)
        for name2 in new_var_names2:
            problem.prob.define(("i", name2), types = lp.VariableType.Binary)

        if exp_obj_name is None:
            exp_obj_name = new_var_id

        bool_ineqs(exp_obj.cont_type(), exp_obj.contain(), new_var_names1, gene_var1, problem)

        bool_ineqs(exp_obj.cont_type(), exp_obj.contain(), new_var_names2, gene_var2, problem)


def bool_ineqs(ctype, containing, names, dict_var, problem):
    '''Input homogenous boolean.Expression (all ANDs or all ORs).  Adds the
    corresponding linear inequalities to the LP problem. No output.

    Args:
        ctype = Identifies container type (AND/OR)
        containing = Opens containers
        names = Names of the variables to replace the logic
        dict_var = dictionary containing already defined expression objects
    '''

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

    if ctype in dict_var.keys():
        Y = dict_var[ctype]
    else:
        Y = 'Y'

    yvar = problem.get_ineq_var(Y)
    ineq_list = []
    for name in names:
        xvar = problem.get_ineq_var(name)
        ineq_list.append(xvar)

    if label == 'or':
        or_group = None
        for ineq_var in ineq_list:
            linear_constraints(problem, 0 <= yvar <= 1)
            linear_constraints(problem, 0 <= ineq_var <= 1)
            #Individual variable inequalities
            orindiv = yvar >= ineq_var
            linear_constraints(problem, orindiv)
            #Container inequalities
            if or_group is None:
                or_group = ineq_var
            else:
                or_group = ineq_var + or_group
        or_cont = yvar <= or_group
        linear_constraints(problem, or_cont)

    if label == 'and':
        and_group = None
        for ineq_var in ineq_list:
            linear_constraints(problem, 0 <= yvar <= 1)
            linear_constraints(problem, 0 <= ineq_var <= 1)
            #Individual variable inequalities
            andindiv = yvar <= ineq_var
            linear_constraints(problem, andindiv)

            #Container inequalities
            if and_group is None:
                and_group = ineq_var
            else:
                and_group = ineq_var + and_group
        and_cont = yvar >= and_group - (N-1)
        linear_constraints(problem, and_cont)


def open_file(self):
    '''Returns the contents of model file in a tuple of dictionaries.
    File Form: tsv format, FOUR Columns: (1) Gene name, (2) Condition 1 Data,
    (3) Condition 2 Data, (4) P-value of the fold change for transition 1->2.'''
    path = self._args.transc_file
    file1 = open(path)
    con1_dict = {}
    con2_dict = {}
    pval_dict = {}

    for row in csv.reader(file1, delimiter=str('\t')):
        try:
            con1_dict[row[0]] = float(row[1])
            con2_dict[row[0]] = float(row[2])
            if float(row[3]) == float(0.0):
                pval_dict[row[0]] = 1e-400
            else:
                pval_dict[row[0]] = float(row[3])
        except:
            print

    return con1_dict, con2_dict, pval_dict


def IDC(dicts):
    '''Used for accessing the list of dictionaries created in open_file()
    Creates a dictionary for the gene ID and a value = [-1, 0, +1] corresponding
    to decreasing, constant, and inreasing expression between the conditions.'''

    con1 = dicts[0]
    con2 = dicts[1]
    pval = dicts[2]
    diff = {}

    for key in con1:
        if con2[key]-con1[key] == 0:
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
            info_list.append(problem.get_flux_var(rxn))
            info[rxn] = info_list
    return info


def add_final_constraints(info_dict, problem, var_ineqvar1, var_ineqvar2):
    '''Takes the output of rxn_info, the LP Problem, and the binary dictionaries
    of each condition.  Adds constraints connecting flux variables, reactions,
    and their flux bounds.'''
    for rxn, info in info_dict.iteritems():

        vmin = info[0]
        vmax = info[1]
        fluxvar = info[2]
        Y = var_ineqvar1[rxn]
        Z = var_ineqvar2[rxn+'.2']

        linear_constraints(problem, fluxvar + (1-Z)*vmax <= vmax)
        linear_constraints(problem, fluxvar + (1-Z)*vmin >= vmin)
