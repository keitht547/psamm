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
from psamm.expression import boolean

from ..command import SolverCommandMixin, MetabolicMixin, Command, CommandError
from .. import fluxanalysis

logger = logging.getLogger(__name__)


class MadeFluxBalance(MetabolicMixin, SolverCommandMixin, Command):
    """Run flux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--loop-removal', help='Select type of loop removal constraints',
            choices=['none', 'tfba', 'l1min'], default='none')
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')

        super(MadeFluxBalance, cls).init_parser(parser)

    def run(self):
        """Run MADE implementation."""
        x = self.parse_dict()
        for key, value in x.iteritems():
            get_root(value)

    def parse_dict(self):
        gene_dict = {}
        for i in self._model.parse_reactions():
            gene_dict[i.id] = i.genes
        return(gene_dict)

def exp_gene_string(G):
    e = boolean.Expression(G)
    #if term in e._root not boolean.Or:
    #print e._root
    #print e.base_tree()
    if type(e.base_tree()) is not boolean.Variable:
        for i in e.base_tree():
            if type(i) is not boolean.Variable:
                for j in i.contain():
                    print(j)
            else:
                print(i)
    else:
        print e.base_tree()

def get_root(G):
    e = boolean.Expression(G)
    #print(e.base_tree())
    if type(e.base_tree()) is boolean.Or:
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
