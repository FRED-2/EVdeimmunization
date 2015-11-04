import numpy
import itertools


class ParetoFilter(object):
    """
        This class implements a pareto filter
        which gets a list of solutions and filters for globally pareto optimal points

        See: Messac et al. 2003
    """
    @staticmethod
    def filter(solutions):
        """
            this function filters a list of solutions for globally pareto optimal points

            @param solutions: A list of solutions
            @type solutions: Solution
            @return: A list of global pareto optimal solutions
        """
        global_sols = []
        m = len(solutions)

        for i in xrange(m):
            is_global = True
            zi = solutions[i]
            for j in xrange(m):
                if i == j:
                    continue

                zj = solutions[j]
                if not zi == zj and all(numpy.greater_equal(zi.objs, zj.objs)):
                    is_global = False
                    break

            if is_global:
                global_sols.append(zi)

        g_sol = []
        for i in xrange(len(global_sols)):
            duplicate = False
            for j in xrange(i+1, len(global_sols)):
                if numpy.allclose(global_sols[i].objs, global_sols[j].objs, rtol=1e-02, atol=1e-06):
                    duplicate = True
                    break

            if not duplicate:
                g_sol.append(global_sols[i])

        return g_sol