
'''Wrapper around Cplex solver'''

import sys
import math
import numbers
from itertools import repeat

import cplex as cp

class CplexSolver(object):
    '''Represents an LP-solver using Cplex'''
    def __init__(self, stream=sys.stderr):
        self._stream = stream

    def create_problem(self):
        '''Create a new LP-problem using the solver'''
        return CplexProblem(stream=self._stream)

class CplexProblem(object):
    '''Represents an LP-problem of a CplexSolver'''

    Maximize = 0
    Minimize = 1

    Continuous = 'C'
    Binary = 'B'
    Integer = 'I'

    def __init__(self, stream=sys.stderr):
        self._cp = cp.Cplex()
        self._cp.set_results_stream(stream)
        self._cp.set_warning_stream(stream)
        self._cp.set_error_stream(stream)
        self._cp.set_log_stream(stream)

        self._variables = set()

    @property
    def cplex(self):
        '''The underlying Cplex object'''
        return self._cp

    def define(self, *names, **kwargs):
        '''Define variable in the problem

        Variables must be defined before they can be accessed by var() or set().
        This function takes keyword arguments lower and upper to define the
        bounds of the variable (default: -inf to inf, where inf is a very large number
        defined by Cplex). The keyword argument types can be used to select the type of
        the variable (Continuous (default), Binary or Integer). Setting any variables
        different than Continuous will turn the problem into an MILP problem.'''

        names = tuple(names)
        lower = kwargs.get('lower', None)
        upper = kwargs.get('upper', None)
        vartype = kwargs.get('types', None)

        # Repeat values if a scalar is given
        if lower is None or isinstance(lower, numbers.Number):
            lower = repeat(lower, len(names))
        if upper is None or isinstance(upper, numbers.Number):
            upper = repeat(upper, len(names))
        if vartype is None or vartype in (CplexProblem.Continuous, CplexProblem.Binary, CplexProblem.Integer):
            vartype = repeat(vartype, len(names))

        # Assign default values
        lower = (-cp.infinity if value is None else value for value in lower)
        upper = (cp.infinity if value is None else value for value in upper)
        vartype = tuple(CplexProblem.Continuous if value is None else value for value in vartype)

        args = { 'names': names, 'lb': tuple(lower), 'ub': tuple(upper) }
        if any(value != CplexProblem.Continuous for value in vartype):
            # Set types only if some are integer (otherwise Cplex will change
            # the solver to MILP).
            args['types'] = vartype

        self._variables.update(names)
        self._cp.variables.add(**args)

    def var(self, name):
        '''Return the variable as an expression'''
        if name not in self._variables:
            raise ValueError('Undefined variable: {}'.format(name))
        return Expression({ name: 1 })

    def set(self, names):
        '''Return the set of variables as an expression

        If any of the variables do not exist in the
        problem, they will be created with the given
        bounds.'''
        names = tuple(names)
        if not self._variables.issuperset(names):
            raise ValueError('Undefined variables: {}'.format(set(names) - self._variables))
        return Expression({ tuple(names): 1 })

    def add_linear_constraints(self, *relations):
        '''Add constraints to the problem

        Each constraint is represented by a Relation, and the
        expression in that relation can be a set expression.'''
        for relation in relations:
            if isinstance(relation, bool):
                # A bool in place of a relation is accepted to mean
                # a relation that does not involve any variables and
                # has therefore been evaluated to a truth-value (e.g
                # '0 == 0' or '2 >= 3').
                if not relation:
                    raise Exception('Unsatisfiable relation added')
            else:
                if relation.sense in (Relation.StrictlyGreater, Relation.StrictlyLess):
                    raise ValueError('Strict relations are invalid in LP-problems: {}'.format(relation))

                expression = relation.expression
                pairs = []
                for value_set in expression.value_sets():
                    ind, val = zip(*((variable, float(value)) for variable, value in value_set))
                    pairs.append(cp.SparsePair(ind=ind, val=val))
                self._cp.linear_constraints.add(lin_expr=pairs, senses=tuple(repeat(relation.sense, len(pairs))),
                                                rhs=tuple(repeat(float(-expression.offset), len(pairs))))

    def set_linear_objective(self, expression):
        '''Set linear objective of problem'''
        self._cp.objective.set_linear(expression.values())

    def set_objective_sense(self, sense):
        '''Set type of problem (maximize or minimize)'''
        if sense == CplexProblem.Minimize:
            self._cp.objective.set_sense(self._cp.objective.sense.minimize)
        elif sense == CplexProblem.Maximize:
            self._cp.objective.set_sense(self._cp.objective.sense.maximize)
        else:
            raise ValueError('Invalid objective sense')

    def solve(self, sense=None):
        '''Solve problem'''
        if sense is not None:
            self.set_objective_sense(sense)
        return self._cp.solve()

    def get_value(self, expression):
        '''Return value of expression'''
        if isinstance(expression, Expression):
            return sum(self._cp.solution.get_values(var)*value for var, value in expression.values())
        elif expression not in self._variables:
            raise ValueError('Unknown expression: {}'.format(expression))
        return self._cp.solution.get_values(expression)

class Expression(object):
    '''Represents a linear expression

    The variables can be ordinary strings which will result
    in a linear expression. If one or more variables are instead
    tuples of strings, then this will be taken to represent a set
    of expressions separately using a different element of the
    tuple.'''
    def __init__(self, variables={}, offset=0):
        self._variables = dict(variables)
        self._offset = offset

    @property
    def offset(self):
        '''Value of the offset'''
        return self._offset

    def variables(self):
        '''Iterator of variables in expression'''
        return self._variables.iterkeys()

    def values(self):
        '''Iterator of variable, value-pairs in expression'''
        return self._variables.iteritems()

    def value_sets(self):
        '''Iterator of expression sets

        This will yield an iterator of variable, value-pairs for
        each expression in the expression set (each equivalent to
        values()). If none of the variables is a set variable then
        a single iterator will be yielded.'''
        count = max(1 if not isinstance(var, tuple) else len(var) for var in self._variables)
        def value_set(n):
            for variable, value in self._variables.iteritems():
                if isinstance(variable, tuple):
                    yield variable[n], value
                else:
                    yield variable, value
        for i in xrange(count):
            yield value_set(i)

    def __add__(self, other):
        '''Add expression with a number or another expression'''
        if isinstance(other, numbers.Number):
            return self.__class__(self._variables, self._offset + other)
        elif isinstance(other, self.__class__):
            result = self.__class__()
            for f in (self._variables, other._variables):
                for var, value in f.iteritems():
                    if var in result._variables:
                        result._variables[var] += value
                    else:
                        result._variables[var] = value
            result._offset = self._offset + other._offset
            return result
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        return self.__class__({ var: value*other for var, value in self._variables.iteritems()}, self._offset*other)

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return self * -1

    def __eq__(self, other):
        '''Return equality relation (equation): self == other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.Equals, self - other)

    def __ge__(self, other):
        '''Return greater-than relation (inequality): self >= other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.Greater, self - other)

    def __le__(self, other):
        '''Return less-than relation (inequality): self <= other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.Less, self - other)

    def __gt__(self, other):
        '''Return strictly greater-than relation (inequality): self > other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.StrictlyGreater, self - other)

    def __lt__(self, other):
        '''Return strictly less-than relation (inequality): self < other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.StrictlyLess, self - other)

    def __str__(self):
        '''Return string representation of expression'''

        def all_terms():
            count_vars = 0
            for name, value in sorted(self._variables.iteritems()):
                if value != 0:
                    count_vars += 1
                    if isinstance(name, tuple):
                        yield '<set>', value
                    else:
                        yield name, value
            if self._offset != 0 or count_vars == 0:
                yield None, self._offset

        terms = []
        for i, spec in enumerate(all_terms()):
            name, value = spec
            if i == 0:
                # First term is special
                if name is None:
                    terms.append('{}'.format(value))
                elif abs(value) == 1:
                    terms.append(name if value > 0 else '-'+name)
                else:
                    terms.append('{}*{}'.format(value, name))
            else:
                prefix = '+' if value >= 0 else '-'
                if name is None:
                    terms.append('{} {}'.format(prefix, abs(value)))
                elif abs(value) == 1:
                    terms.append('{} {}'.format(prefix, name))
                else:
                    terms.append('{} {}*{}'.format(prefix, abs(value), name))
        return ' '.join(terms)

    def __repr__(self):
        return '<Expression \'{}\'>'.format(str(self))

class Relation(object):
    '''Represents a binary relation (equation or inequality)

    Relations can be equalities or inequalities. All relations
    of this type can be represented as a left-hand side expression
    and the type of relation. In this representation, the right-hand
    side is always zero.'''

    Equals = 'E'
    Greater = 'G'
    Less = 'L'
    StrictlyGreater = 'SG'
    StrictlyLess = 'SL'

    SYMBOL = {
        Equals: '==',
        Greater: '>=',
        Less: '<=',
        StrictlyGreater: '>',
        StrictlyLess: '<'
    }

    def __init__(self, sense, expression):
        self._sense = sense
        self._expression = expression

    @property
    def sense(self):
        '''Type of relation (equality or inequality)

        Can be one of Equal, Greater or Less, or one of the
        strict relations, StrictlyGreater or StrictlyLess.'''
        return self._sense

    @property
    def expression(self):
        '''Left-hand side expression'''
        return self._expression

    def __str__(self):
        '''Convert relation to string representation'''
        return '{} {} 0'.format(str(self._expression), Relation.SYMBOL[self._sense])

    def __repr__(self):
        return '<Relation \'{}\'>'.format(str(self))


def convex_cardinality_relaxed(f, epsilon=1e-5):
    '''Transform L1-norm optimization function into approximate cardinality optimization

    The given function must optimize a convex problem with
    a weighted L1-norm as the objective. The transformed function
    will apply the iterated weighted L1 heuristic to approximately
    optimize the cardinality of the solution. This method is
    described by S. Boyd, "L1-norm norm methods for convex cardinality
    problems." Lecture Notes for EE364b, Stanford University, 2007.
    Available online at www.stanford.edu/class/ee364b/.

    The given function must take an optional keyword parameter weights
    (dictionary), and the weights must be set to one if not specified.
    The function must return the non-weighted solution as an iterator
    over (identifier, value)-tuples, either directly or as the first
    element of a tuple.'''

    def convex_cardinality_wrapper(*args, **kwargs):
        def dict_result(r):
            if isinstance(r, tuple):
                return dict(r[0])
            return dict(r)

        # Initial run with default weights
        full_result = f(*args, **kwargs)
        result = dict_result(full_result)

        def update_weight(value):
            return 1/(epsilon + abs(value))

        # Iterate until the difference from one iteration to
        # the next is less than epsilon.
        while True:
            weights = { identifier: update_weight(value) for identifier, value in result.iteritems() }
            kwargs['weights'] = weights

            last_result = result
            full_result = f(*args, **kwargs)
            result = dict_result(full_result)

            delta = math.sqrt(sum(pow(value - last_result[identifier], 2) for identifier, value in result.iteritems()))
            if delta < epsilon:
                break

        if isinstance(full_result, tuple):
            return (result.iteritems(),) + full_result[1:]
        return result.iteritems()

    return convex_cardinality_wrapper