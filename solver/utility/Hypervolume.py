import heapq
import Polygon


class HyperVolume(object):
    """
        this class can be used to measure the quality of a pareto approximation
    """

    @staticmethod
    def calc_hypervol(solutions, init_rec):
        """
            calculates the hypervolume indicator according to Zitler and Thiele


            @param init_rec: the rectangle spanning the search space (represented by two diagonal vertices)
            @return: the hypervolume and the Polygon spanned by the pareto point and the search space

            @assertion: self._solution is encoded as heapq
        """
        sorted_sol = [s.objs for s in heapq.nsmallest(len(solutions), solutions)]
        first = sorted_sol[0]
        last = sorted_sol[-1]

        #add now intersection points between solution poly and init_rect
        sorted_sol.append([first[0], init_rec[0][1]])
        sorted_sol.append([init_rec[1][0], last[1]])
        sorted_sol.append([init_rec[1][0], init_rec[0][1]])

        p = Polygon.Polygon(sorted_sol)
        return p.area(), p

    @staticmethod
    def calc_adj_hypervol(hyper_poly, solutions, sol_rec):
        """
            calculates the adjusted hypervolume according to Boland et al. 2013

            @param hyper_poly: the polygon spanned by the pareto points and the search space
            @param sol_rec: A dictionary holding rectangles for which it is proven that they don't contain
            unrecognized nondominating points
            @return: The calculated area and the spanning polygon

            @assertion: self._solution is encoded as heapq!!!
        """
        sorted_sol = [s.objs for s in heapq.nsmallest(len(solutions), solutions)]
        for i in xrange(len(sorted_sol)-1):
            z_t = sorted_sol[i]
            z_b = sorted_sol[i+1]
            #if one could prove that the spanning rectangle does not contain further nondominating points
            #exclude from area calculation
            if (z_t, z_b) in sol_rec:
                continue
            #now add the spanning rectangle to the hyper_poly
            hyper_poly += Polygon.Polygon([z_t, [z_b[0], z_t[1]], z_b, [z_t[0], z_b[1]]])
        return hyper_poly.area(), hyper_poly

    @staticmethod
    def calc_hypervol_gap(solutions, init_rec, empty_rec):
        """
            calculates the percentage of the gap between hypervolume and adjusted hypervolume
        """
        hyper_vol, hyper_poly = HyperVolume.calc_hypervol(solutions, init_rec)
        adj_hyper_vol, ady_hyper_poly = HyperVolume.calc_adj_hypervol(hyper_poly, solutions, empty_rec)
        return (adj_hyper_vol - hyper_vol) / hyper_vol