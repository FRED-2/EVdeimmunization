from __future__ import division
import cplex
import argparse
import numpy
import heapq
import subprocess
import time
import os
import socket
import multiprocessing as mp
import itertools as itr
from multiprocessing.managers import SyncManager
from utility import HyperVolume
from utility import ParetoFilter




class RectangleEpsilonGridManager(object):
    """
        This is prarallel implementation of
        the epsilon contraint method with precomputed boundaries

        CPLEX floating point accuracy is b default ~1e-6
    """

    def __init__(self, port=6881, authkey="rectangle", constraints=None, output=None, verbose=0):


        self.solutions = []
        self.biob_cons = ["z2_cons", "z1_cons"] if constraints is None else constraints
        self._hypervol = HyperVolume()
        self._manager = self.__make_manager_server(port, authkey)
        #concurrent stuff
        self.task_q = self._manager.get_task_q()
        self.done_q = self._manager.get_done_q()
        self.empty_rectangles = []
        self.output = output
        self.verbose = verbose

    @staticmethod
    def __isclose(a, b, rel_tol=1e-09, abs_tol=0.00001):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    def __make_manager_server(self, port, authkey):
        """
         Starts a manager server
        :return: Manager object
        """

        job_q = mp.JoinableQueue()
        result_q = mp.Queue()
        result_eps_q = mp.Queue()
        task_q = mp.JoinableQueue()

        # This is based on the examples in the official docs of multiprocessing.
        # get_{job|result}_q return synchronized proxies for the actual Queue
        # objects.
        class JobQueueManager(mp.managers.SyncManager):
            pass

        JobQueueManager.register('get_task_q', callable=lambda: job_q)
        JobQueueManager.register('get_done_q', callable=lambda: result_q)
        JobQueueManager.register('get_eps_done_q', callable=lambda: result_eps_q)
        JobQueueManager.register('get_task_q', callable=lambda: task_q)
        manager = JobQueueManager(address=('', port), authkey=authkey)

        manager.start()
        print 'Server started at %s:%s. Authentication: %s' % (socket.gethostbyname(socket.gethostname()),port,authkey)
        return manager

    def solve_grid(self, nof_sol=1):


        solutions = []
        #init problems to solve
        self.task_q.put((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
        self.task_q.put((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

        self.task_q.join()

        b = [None, None]
        warm = [None, None]
        while not self.done_q.empty():
            pos, sol, warmstart, origin_rect = self.done_q.get()
            print "Edge Point ", sol.objs
            solutions.append(sol)
            b[pos] = sol.objs
            warm[pos] = warmstart[pos]

        #temporal save
        #pcl.dump(solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

        delta = 1/(nof_sol)
        diff = b[1][1] - b[0][1]
        alphas = [delta*i*diff for i in xrange(1, nof_sol+1)]

        c = 0
        for i in xrange(nof_sol):
            print "Bounds: ", b[0][1]+alphas[i]
            self.task_q.put((0, 1, b[0][1]+alphas[i], warm, tuple(b)))
            c += 1

        while c:
            _, sol, _, _ = self.done_q.get()
            print sol
            c -= 1
            solutions.append(sol)
            #pcl.dump(solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

        print "solved all epsilon grids"
        return solutions


    def solve_rectangle(self, init_recs=None, rel_tol=1e-05, abs_tol=0.0001):
        """
            solve rectangle spliting after intial deterimnation fo rectangles
        """
        eps=abs_tol
        isclose =RectangleEpsilonGridManager.__isclose
        task_count = 0
        if init_recs:
            b = [init_recs[0].objs, init_recs[-1].objs]
            self.solutions.extend(init_recs)
            for i in xrange(len(init_recs)-1):
                task_count += 1
                zi = init_recs[i]
                zj = init_recs[i+1]
                warm = [zi.warm_start, zj.warm_start]
                rec = [zi.objs, zj.objs]
                rec_b = (rec[0][1]+rec[1][1])/2
                self.task_q.put((0, 1, rec_b, warm, tuple(rec)))
        else:
            task_count = 2
            #init problems to solve
            self.task_q.put((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
            self.task_q.put((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

            self.task_q.join()


            b = [None, None]
            warm = [None, None]
            while task_count:
                pos, sol, warmstart, origin_rect = self.done_q.get()
                if pos is None:
                    continue
                task_count -= 1
                self.solutions.append(sol)
                #pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                b[pos] = sol.objs
                warm[pos] = warmstart[pos]

            rec_b = (b[0][1]+b[1][1])/2
            task_count = 1
            self.task_q.put((0, 1, rec_b, warm, b))

        while task_count:

            pos, sol, warm, origin_rect = self.done_q.get()

            #lexmin2
            if pos:
                if not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[0])): 
                   if not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[1])):
                    #if not numpy.allclose(sol.objs, origin_rect[0], rtol=1e-04, atol=1e-02):
                        rec = (origin_rect[0], sol.objs)
                        self.solutions.append(sol)
                        
                        if abs(rec[0][0]-rec[1][0]) > eps:
                            rec_b = (rec[0][1]+rec[1][1])/2
                            self.task_q.put((0, 1, rec_b, warm, rec))
                            task_count += 1
            #lexmin1
            else:
                rec_t = sol.objs[0]-eps
                self.task_q.put((1, 0, rec_t, warm, origin_rect))
                task_count += 1
                if not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[1])):
                   if not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[0])):
                    #if not numpy.allclose(sol.objs, origin_rect[1], rtol=1e-04, atol=1e-02):
                        rec = (sol.objs, origin_rect[1])
                        self.solutions.append(sol)

                        if abs(rec[0][0]-rec[1][0]) > eps:
                            rec_b = (rec[0][1]+rec[1][1])/2
                            self.task_q.put((0, 1, rec_b, warm, rec))
                            task_count += 1

            task_count -= 1

            if self.verbose:
                print "Current Rectangle ", origin_rect
                print "Solution ", sol
                print "Lexmin:",pos
                print "Solutions ", len(self.solutions)
                print "Rectangle Size: ", abs(origin_rect[0][0]-origin_rect[1][0])
                print "Tasks still running: ", task_count

        return self.solutions


    def approximate_rectangle(self, gap, init_recs=None, rel_tol=1e-05, abs_tol=0.0001):
        eps = abs_tol
        verbose = self.verbose
        isclose = RectangleEpsilonGridManager.__isclose
        task_count = 0
        proof_not_empty_rec = {}
        empty_rect = set()

        if init_recs:
            init_rect = [init_recs[0].objs, init_recs[-1].objs]
            self.solutions.extend(init_recs)
            for i in xrange(len(init_recs)-1):
                task_count += 1
                zi = init_recs[i]
                zj = init_recs[i+1]
                warm = [zi.warm_start, zj.warm_start]
                rec = [zi.objs, zj.objs]
                rec_b = (rec[0][1]+rec[1][1])/2
                self.task_q.put((0, 1, rec_b, warm, tuple(rec)))
        else:
            task_count += 2
            self.task_q.put((0, 1, cplex.infinity, [None, None], ((None, None), (None, None))))
            self.task_q.put((1, 0, cplex.infinity, [None, None], ((None, None), (None, None))))

            self.task_q.join()

            init_rect = [None, None]
            warm = [None, None]
            while task_count:
                pos, sol, warmstart, origin_rect = self.done_q.get()

                if pos is None:
                    continue
                task_count -= 1
                heapq.heappush(self.solutions, sol)
                #pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                init_rect[pos] = sol.objs
                warm[pos] = warmstart[pos]
            init_rect = tuple(init_rect)
            rec_b = (init_rect[0][1]+init_rect[1][1])/2
            task_count += 1
            self.task_q.put((0, 1, rec_b, warm, init_rect))

        while task_count > 0:
            pos, sol, warm, origin_rect = self.done_q.get()

            if verbose:
                print
                print "Current Rectangle ", origin_rect
                print "Solution: ", sol
                print "Tasks still running: ", task_count

            if pos:
                if not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[1])) and not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[0])):
                    rec = (origin_rect[0], sol.objs)
                    heapq.heappush(self.solutions, sol)
                    #pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))
                    if abs(rec[0][0]-rec[1][0]) > eps:
                        rec_b = (rec[0][1]+rec[1][1])/2
                        self.task_q.put((0, 1, rec_b, warm, rec))
                        task_count += 1

                    if origin_rect in proof_not_empty_rec:
                        print "Rectanlge ",(sol.objs, proof_not_empty_rec[origin_rect])," is empty"
                        empty_rect.add((sol.objs, proof_not_empty_rec[origin_rect]))
                        del proof_not_empty_rec[origin_rect]

                        cur_gap = self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
                        print "Current Hypervol. gap is: ", cur_gap
                        if cur_gap < gap:
                           break

                    else:
                        #empty_rect.add((sol.objs, origin_rect[1]))
                        proof_not_empty_rec[origin_rect] = sol.objs

                else:

                    if origin_rect in proof_not_empty_rec:
                        print "Rectanlge ",(origin_rect[0], proof_not_empty_rec[origin_rect])," is empty"
                        empty_rect.add((origin_rect[0], proof_not_empty_rec[origin_rect]))
                        del proof_not_empty_rec[origin_rect]

                        cur_gap = self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
                        print "Current Hypervol. gap is: ", cur_gap
                        if cur_gap < gap:
                           break
                    else:
                        proof_not_empty_rec[origin_rect] = origin_rect[0]

            #lexmin1
            else:
                rec_t = sol.objs[0]-eps
                self.task_q.put((1, 0, rec_t, warm, origin_rect))
                task_count += 1

                if not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[1])) and not all(isclose(a,b, rel_tol=rel_tol, abs_tol=abs_tol) for a,b in itr.izip(sol.objs, origin_rect[0])):
                    rec = (sol.objs, origin_rect[1])
                    heapq.heappush(self.solutions, sol)

                    if abs(rec[0][0]-rec[1][0]) > eps:
                        rec_b = (rec[0][1]+rec[1][1])/2
                        self.task_q.put((0, 1, rec_b, warm, rec))
                        task_count += 1

                    if origin_rect in proof_not_empty_rec:
                        empty_rect.add((proof_not_empty_rec[origin_rect], sol.objs))
                        del proof_not_empty_rec[origin_rect]

                        cur_gap = self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
                        if verbose:
                            print "Rectanlge ",(proof_not_empty_rec[origin_rect], origin_rect[1])," is empty"
                            print "Current Hypervol. gap is: ", cur_gap
                        if cur_gap < gap:
                           break
                    else:
                        proof_not_empty_rec[origin_rect] = sol.objs
                else:
                    if origin_rect in proof_not_empty_rec:
                        empty_rect.add((proof_not_empty_rec[origin_rect], origin_rect[1]))
                        del proof_not_empty_rec[origin_rect]

                        cur_gap = self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
                        if verbose:
                            print "Rectanlge ",(proof_not_empty_rec[origin_rect], sol.objs)," is empty"
                            print "Current Hypervol. gap is: ", cur_gap
                        if cur_gap < gap:
                           break

                    else:
                        proof_not_empty_rec[origin_rect] = origin_rect[1]

            task_count -= 1


        if verbose:
            print "Last hypervol gap: ", self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)

        return self.solutions


    def solve(self, nof_sol=3, rel_tol=1e-05, abs_tol=0.0001):
        sols = self.solve_grid(nof_sol=nof_sol)

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols,rel_tol=rel_tol, abs_tol=abs_tol))

        if self.verbose:
            print "Filtered pareto points"
            for s in sols:
                print s

        #pcl.dump(sols, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))
        self.solve_rectangle(init_recs=sols,rel_tol=rel_tol, abs_tol=abs_tol)
        for _ in xrange(nof_sol):
            self.task_q.put(("DONE",None,None,None,None))
        self._manager.shutdown()
        return self.solutions

    def approximate(self, gap, nof_sol=3, rel_tol=1e-05, abs_tol=0.0001):
        sols = self.solve_grid(nof_sol=nof_sol)

        #filter non pareto points
        sols = ParetoFilter.filter(sorted(sols),rel_tol=rel_tol, abs_tol=abs_tol)

        if self.verbose:
            print "Filtered pareto points"
            for s in sols:
                print s

            print "\n","Starting rectangle refinement\n"

        curr_gap = HyperVolume.calc_hypervol_gap(sols, [sols[0].objs, sols[-1].objs], [])
        if curr_gap < gap:
            return sols

        self.approximate_rectangle(gap, init_recs=sols, rel_tol=1e-05, abs_tol=0.0001)
        for _ in xrange(nof_sol):
            self.task_q.put(("DONE",None,None,None,None))

        self._manager.shutdown()
        return self.solutions


if __name__ == "__main__":
    import cPickle as pcl


    parser = argparse.ArgumentParser(description=' Rectangle Manager implementation')
    parser.add_argument('--worker','-w',
                      required=False,
                      type=int,
                      default=3,
                      help="Number number of worker processes to be started (default:3)")
    parser.add_argument('--port','-p',
                      required=True,
                      type=int,
                      help="Port")
    parser.add_argument('--approximate','-a',
                      type=float,
                      default=0.009,
                      help="Bound on approximation (default:0.009)")
    parser.add_argument('--key','-k',
                      type=str,
                      default="rectangle",
                      help="Authentication key (default: rectangle)")
    parser.add_argument('--relTol','-rel',
                      default=0.00001,
                      type=float,
                      required=False,
                      help="The relative tollerance for floating point comparison")
    parser.add_argument('--absTol','-abs',
                      default=0.001,
                      type=float,
                      required=False,
                      help="The absolut tollerance for floating point comparison, also used as epsilon in the model")
    parser.add_argument('--output','-o',
                      required=True,
                      help="Solution output as pickel")
    parser.add_argument('--resolve','-r',
                      required=False,
                      help="Reinitialize with partial solution")
    parser.add_argument('--verbose','-v',
                      default=0,
                      type=int,
                      required=False,
                      help="Verbosity (default 0)")

    args = parser.parse_args()
    manager = RectangleEpsilonGridManager(port=args.port, authkey=args.key, output=args.output, verbose=args.verbose)
    if args.resolve:
        f = open(args.resolve, "r")
        sol = sorted(pcl.load(f))
        s = time.time()
        sols=manager.approximate_rectangle(args.approximate, init_recs=sol, nof_sol=args.worker, 
                                           rel_tol=args.relTol, abs_tol=args.absTol)   
        manager._manager.shutdown() 
        e=time.time()
    else:    
        s = time.time()
        sols = manager.solve(nof_sol=args.worker,rel_tol=args.relTol, abs_tol=args.absTol)
        #sols = manager.approximate(args.approximate, nof_sol=args.worker)
        #sols = manager.solve_rectangle(rel_tol=args.relTol, abs_tol=args.absTol)
        e= time.time()
        filtered_sols = ParetoFilter.filter(sols)
    print "Solution time ", e-s
    for s in sols:
        print s
    print "NOF: ",len(sols), len(filtered_sols)
    #pcl.dump(sols, open(args.output, "w"))