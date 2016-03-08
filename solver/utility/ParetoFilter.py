import numpy
import itertools


class ParetoFilter(object):
    """
        This class implements a pareto filter
        which gets a list of solutions and filters for globally pareto optimal points

        See: Messac et al. 2003
    """
    @staticmethod
    def filter(solutions, rel_tol=1e-05, abs_tol=0.0001):
        """
            this function filters a list of solutions for globally pareto optimal points

            @param solutions: A list of solutions
            @type solutions: Solution
            @return: A list of global pareto optimal solutions
        """
        solutions.sort()
        isclose=lambda a,b: abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
        g_sol = []
        for i in xrange(len(solutions)-1):
            duplicate = False
            for j in xrange(i+1, len(solutions)):
                if all(isclose(a,b) for a,b in itertools.izip(solutions[i].objs,solutions[j].objs)):
                    print "z1=",solutions[i].objs," equal to ",solutions[j].objs
                    duplicate = True
                    break

            if not duplicate:
                g_sol.append(solutions[i])
        g_sol.append(solutions[-1])

        m = len(g_sol)
        final=[]
        for i in xrange(m):
            is_global = True
            zi = g_sol[i]
            for j in xrange(m):
                if i == j:
                    continue

                zj = g_sol[j]
                if all(numpy.greater_equal(zi.objs, zj.objs)):
                    print "z1=",zi.objs," is exlcuded due to:",zj.objs
                    is_global = False
                    break

            if is_global:
                final.append(zi)

        return final