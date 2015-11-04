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
from multiprocessing.managers import SyncManager
from utility import HyperVolume
from utility import ParetoFilter




class RectangleEpsilonGridManager(object):
    """
        This is prarallel implementation of
        the epsilon contraint method with precomputed boundaries
    """
    EPS = 1e-3

    def __init__(self, port=6881, authkey="rectangle", constraints=None, output=None):


        self.solutions = []
        self.biob_cons = ["z2_cons", "z1_cons"] if constraints is None else constraints
        self._hypervol = HyperVolume()
        self._manager = self.__make_manager_server(port, authkey)
        #concurrent stuff
        self.task_q = self._manager.get_task_q()
        self.done_q = self._manager.get_done_q()
        self.empty_rectangles = []
        self.output = output


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
        pcl.dump(solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

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
            pcl.dump(solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

        print "solved all epsilon grids"
        return solutions


    def solve_rectangle(self, init_recs=None):
        """
            solve rectangle spliting after intial deterimnation fo rectangles
        """

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
                rec_b = 0.5*(rec[0][1]+rec[1][1])
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
                pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                b[pos] = sol.objs
                warm[pos] = warmstart[pos]

            rec_b = 0.5*(b[0][1]+b[1][1])
            task_count = 1
            self.task_q.put((0, 1, rec_b, warm, b))

        while task_count:

            pos, sol, warm, origin_rect = self.done_q.get()
            print "Current Rectangle ", origin_rect
            print "Solution ", sol
            print "Solutions ", self.solutions

            #lexmin2
            if pos:
                if not numpy.allclose(sol.objs, origin_rect[0], rtol=1e-03, atol=1e-04):
                    rec = (origin_rect[0], sol.objs)
                    self.solutions.append(sol)

                    pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    self.task_q.put((0, 1, rec_b, warm, rec))
                    task_count += 1
            #lexmin1
            else:
                rec_t = sol.objs[0]-RectangleEpsilonGridManager.EPS
                self.task_q.put((1, 0, rec_t, warm, origin_rect))
                task_count += 1
                if not numpy.allclose(sol.objs, origin_rect[1], rtol=1e-03, atol=1e-04):
                    rec = (sol.objs, origin_rect[1])
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    self.solutions.append(sol)

                    pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                    self.task_q.put((0, 1, rec_b, warm, rec))
                    task_count += 1

            task_count -= 1
            print "Tasks still running: ", task_count



    def approximate_rectangle(self, gap, init_recs=None):
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
                rec_b = 0.5*(rec[0][1]+rec[1][1])
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
                pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                init_rect[pos] = sol.objs
                warm[pos] = warmstart[pos]
            init_rect = tuple(init_rect)
            rec_b = 0.5*(init_rect[0][1]+init_rect[1][1])
            task_count += 1
            self.task_q.put((0, 1, rec_b, warm, init_rect))

        while task_count > 0:

            pos, sol, warm, origin_rect = self.done_q.get()
            print
            print
            print "Current Rectangle ", origin_rect
            print "Solution: ", sol

            if pos:
                if not numpy.allclose(sol.objs, origin_rect[0], rtol=1e-03, atol=1e-04):
                    rec = (origin_rect[0], sol.objs)
                    heapq.heappush(self.solutions, sol)
                    pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    self.task_q.put((0, 1, rec_b, warm, rec))
                    task_count += 1

                    if origin_rect in proof_not_empty_rec:
                        print "Rectanlge ",(sol.objs, proof_not_empty_rec[origin_rect])," is empty"
                        empty_rect.add((sol.objs, proof_not_empty_rec[origin_rect]))
                        del proof_not_empty_rec[origin_rect]
                    else:
                        #empty_rect.add((sol.objs, origin_rect[1]))
                        proof_not_empty_rec[origin_rect] = sol.objs

                else:

                    if origin_rect in proof_not_empty_rec:
                        print "Rectanlge ",(origin_rect[0], proof_not_empty_rec[origin_rect])," is empty"
                        empty_rect.add((origin_rect[0], proof_not_empty_rec[origin_rect]))
                        del proof_not_empty_rec[origin_rect]
                    else:
                        proof_not_empty_rec[origin_rect] = origin_rect[0]

            #lexmin1
            else:
                rec_t = sol.objs[0]-RectangleEpsilonGridManager.EPS
                self.task_q.put((1, 0, rec_t, warm, origin_rect))
                task_count += 1

                if not numpy.allclose(sol.objs, origin_rect[1], rtol=1e-03, atol=1e-04):
                    rec = (sol.objs, origin_rect[1])
                    rec_b = 0.5*(rec[0][1]+rec[1][1])
                    heapq.heappush(self.solutions, sol)
                    pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

                    self.task_q.put((0, 1, rec_b, warm, rec))
                    task_count += 1
                    if origin_rect in proof_not_empty_rec:
                        print "Rectanlge ",(proof_not_empty_rec[origin_rect], sol.objs)," is empty"
                        empty_rect.add((proof_not_empty_rec[origin_rect], sol.objs))
                        del proof_not_empty_rec[origin_rect]

                    else:
                        proof_not_empty_rec[origin_rect] = sol.objs
                else:
                    if origin_rect in proof_not_empty_rec:
                        print "Rectanlge ",(proof_not_empty_rec[origin_rect], origin_rect[1])," is empty"
                        empty_rect.add((proof_not_empty_rec[origin_rect], origin_rect[1]))
                        del proof_not_empty_rec[origin_rect]

                    else:
                        proof_not_empty_rec[origin_rect] = origin_rect[1]

            task_count -= 1
            print "Tasks still running: ", task_count

            #if desired approximation quality is reached terminate and return results
            cur_gap = self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
            print "Current Hypervol. gap is: ", cur_gap

            if numpy.allclose(gap, cur_gap, rtol=1e-03, atol=1e-04) or cur_gap < gap:
                break
        print
        print
        print "Last hypervol gap: ", self._hypervol.calc_hypervol_gap(self.solutions, init_rect, empty_rect)
        pcl.dump(self.solutions, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

        return self.solutions


    def solve(self, nof_sol=3):
        sols = self.solve_grid(nof_sol=nof_sol)

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols))

        print "Filtered pareto points"
        for s in sols:
            print s

        pcl.dump(sols, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

        self.solve_rectangle(init_recs=sols)
        for _ in xrange(nof_sol):
            self.task_q.put(("DONE",None,None,None,None))
        self._manager.shutdown()
        return self.solutions

    def approximate(self, gap, nof_sol=3):
        sols = self.solve_grid(nof_sol=nof_sol)

        #filter non pareto points
        sols = sorted(ParetoFilter.filter(sols))

        print "\nFiltered pareto points"
        for s in sols:
            print s

        print "\n","Starting rectangle refinement\n"
        pcl.dump(sols, open(".".join(args.output.split(".")[:-1])+"_temp.pcl", "w"))

        curr_gap = HyperVolume.calc_hypervol_gap(sols, [sols[0].objs, sols[-1].objs], [])
        if numpy.allclose(gap, curr_gap, rtol=1e-03, atol=1e-04) or curr_gap < gap:
            return sols

        self.approximate_rectangle(gap, init_recs=sols)
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
    parser.add_argument('--output','-o',
                      required=True,
                      help="Solution output as pickel")
    parser.add_argument('--resolve','-r',
                      required=False,
                      help="Reinitialize with partial solution")

    args = parser.parse_args()
    manager = RectangleEpsilonGridManager(port=args.port, authkey=args.key, output=args.output)
    if args.resolve:
        f = open(args.resolve, "r")
        sol = sorted(pcl.load(f))
        print sol
        s = time.time()
        sols=manager.approximate_rectangle(args.approximate, init_recs=sol, nof_sol=args.worker)   
        manager._manager.shutdown() 
        e=time.time()
    else:    
        s = time.time()
        sols = manager.approximate(args.approximate, nof_sol=args.worker)
        e= time.time()

    print "Solution time ", e-s
    for s in sols:
        print s
    pcl.dump(sols, open(args.output, "w"))