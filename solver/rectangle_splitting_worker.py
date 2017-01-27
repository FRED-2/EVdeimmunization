from __future__ import division
import itertools
import sys
import cplex
import numpy
import argparse
import time 
from multiprocessing.managers import SyncManager
from utility import Solution

import time

class RectangleSplittingWorker(object):

    def __init__(self, z1_name, z2_name, biob_cons, inter_vars, port, authkey, ip, nof_cpu=6, has_constraints=False, verbose=0):
        z1 = cplex.Cplex(z1_name)
        z2 = cplex.Cplex(z2_name)
        z1.parameters.threads.set(int(nof_cpu))
        z2.parameters.threads.set(int(nof_cpu))

        if not verbose:
            z1.set_log_stream(None)
            z1.set_error_stream(None)
            z1.set_warning_stream(None)
            z1.set_results_stream(None)
            z2.set_log_stream(None)
            z2.set_error_stream(None)
            z2.set_warning_stream(None)
            z2.set_results_stream(None)

        self.verbose=verbose
        self.manager = self.__make_client_manager(port, authkey, ip)
        self._models = (z1, z2)
        self._changeable_constraints = biob_cons
        self._variables = (z1.variables.get_names(),z2.variables.get_names())
        self._inter_variables = filter(lambda x:  x[0] in inter_vars, self._variables[0])
        self.task_q = self.manager.get_task_q()
        self.done_q = self.manager.get_done_q()

        if verbose:
            print "Modify model for solving"
        #c=1
        
        '''
         takes longer than external construction of constraints due to internal string matching from one
         single-objective problem to the other
        '''
        s = time.time()
        if not has_constraints:
           
            allclose = numpy.allclose
            z1_obj_val = { k:v for k,v in itertools.izip(self._variables[0],z1.objective.get_linear()) if not allclose(v, 0)}
            z2_obj_val = { k:v for k,v in itertools.izip(self._variables[1],z2.objective.get_linear()) if not allclose(v, 0)}    


            z1.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z2_obj_val.keys(), val=z2_obj_val.values())], senses=["L"], rhs=[0.0], range_values=[0], names=[biob_cons[0]])
            z2.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=z1_obj_val.keys(), val=z1_obj_val.values())], senses=["L"], rhs=[0.0], names=[biob_cons[1]])

        if verbose:
            e = time.time()
            print "Construction took: ",e-s
        #sys.exit()
        #run the worker
        self.run()

    def __make_client_manager(self, port, authkey, ip):
        """
        generates manager and connects to master manager
        :param port:
        :param authkey:
        :return: manager
        """
        class ServerQueueManager(SyncManager):
            pass

        ServerQueueManager.register('get_task_q')
        ServerQueueManager.register('get_done_q')

        manager = ServerQueueManager(address=(ip, port), authkey=authkey)
        manager.connect()

        if self.verbose:
            print 'Client connected to %s:%s' % (ip, port)
        return manager

    def _lexmin(self, z1_idx, z2_idx, boundary,  warmstart=None, effort_level=0):
        """
            lexicographic optimization based on input values

            :param z1_idx, z2_idx: integers defining which of the two objectives is optimized first
            :param z1_bound, boundary: floating point of new boundary of the
            :param warmstart (optional): is an optional feature to use a warmstart for z1
            :return Solution, values_of_all_variables_for_warm_start


            :QUESTION:Dont know if one should first delete old MIP starts or just keep adding them
        """

        z1 = self._models[z1_idx]
        
        if warmstart:
            z1.MIP_starts.delete() #this is questionable
            z1.MIP_starts.add([self._variables[z1_idx], warmstart], effort_level)
        z1.linear_constraints.set_rhs(self._changeable_constraints[z1_idx], boundary)
        z1.solve()
        z1_hat = z1.solution.get_objective_value()
        z1_hat_values = z1.solution.get_values()

        #lexicographical solution second model
        z2 = self._models[z2_idx]

        z2.linear_constraints.set_rhs(self._changeable_constraints[z2_idx], z1_hat)
        #set warm start from first objective
        z2.MIP_starts.add([self._variables[z1_idx], z1_hat_values], z2.MIP_starts.effort_level.auto)
        z2.solve()
        z2_hat = z2.solution.get_objective_value()
        z2_hat_values = z2.solution.get_values()

        #generate solution object
        objs = [0, 0]
        objs[z1_idx] = z1_hat
        objs[z2_idx] = z2_hat

        inter_vars = {}
        if len(self._inter_variables) > 0:
            inter_vars = {k:v for v, k in itertools.izip(z2.solution.get_values(self._inter_variables),
                                                         self._inter_variables)
                          if 0.9<=v<=1.02 }
        s = Solution(objs, inter_vars)

        return s, z2_hat_values

    def run(self):
        while True:
            z1_idx, z2_idx, boundary, warmst, rectangle = self.task_q.get()

            if z1_idx == "DONE":
                self.task_q.task_done()
                sys.exit()

            if self.verbose:
                print "solving ",z1_idx, " in recangle ", rectangle, " with boundary ",boundary

            try:
                sol, warm = self._lexmin(z1_idx, z2_idx, boundary,  warmstart=warmst[0], effort_level=0)
                if z1_idx:
                    self.done_q.put((z1_idx, sol, [warmst[0], warm], rectangle))
                else:
                    self.done_q.put((z1_idx, sol, [warm, warmst[1]], rectangle))
            except Exception as e:
                print "An error has occured or the server is shut-down: ",e
                print "Shutting down...."
                sys.exit()            
            self.task_q.task_done()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' Rectangle Worker Grid implementation')
    parser.add_argument('--input','-i',
                      required=True,
                      nargs=2,
                      help="model files ")

    parser.add_argument('--masterip','-m',
                      required=True,
                      help="The IP of the master node"
                      )
    parser.add_argument('--port','-p',
                      type=int,
                      required=True,
                      help="port to connect"
                      )
    parser.add_argument('--authkey','-a',
                      default="rectangle",
                      required=False,
                      help="authentication key"
                      )
    parser.add_argument('--verbose','-v',
                      default=0,
                      type=int,
                      required=False,
                      help="verbosity"
                      )
    parser.add_argument('--threads','-t',
                      type=int,
                      required=True,
                      help="nof of core"
                      )

    args = parser.parse_args()
    worker = RectangleSplittingWorker(args.input[0], args.input[1],
                                      ["z2_cons", "z1_cons"],["x","y"], args.port, args.authkey,
                                      args.masterip, args.threads, False, args.verbose)

    sys.exit()

