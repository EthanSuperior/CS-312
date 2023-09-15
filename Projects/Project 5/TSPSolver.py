#!/usr/bin/python3
import copy
import random
from turtle import clone

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy( self,time_allowance=60.0 ):
		"""
		This solves the TSP problem by using a greedy approach with the Nearest Neighbor (NN) algorithm.
		It chooses a random city, and then finds the next city it has not visited by getting the cheapest cost.
		The normal approach to this only runs on 1 index, but mine searches up to all the nodes to find a solution
		Complexity:
			Time: O(n^3)
			Space: O(n)

		:param time_allowance: Maximum time allowed to generate a solution for.
		:return: A TSP path it one was found otherwise np.inf
		"""
		'''Initialize the cites from the scenario'''
		results, cities = {}, self._scenario.getCities()
		'''Initialize some of the data we need for results'''
		ncities, count, bssf = len(cities), 0, None
		start_time = time.time()
		''' Random Permutation for the possible starting indexes;
		 	this allows multiple NN searches to be run with no overlaps
		 		Time: O(n) Space: O(n)
		 	'''
		for k in np.random.permutation(ncities):
			'''Terminate if we are taking too long'''
			if time.time()-start_time >= time_allowance: break
			'''Create a list of city indexes to search through'''
			unvisited_cities = list(range(ncities))
			'''Initialize the starting route to random starting index and the from_city to the respective city'''
			from_city, route = cities[k], [unvisited_cities.pop(k)]
			'''Search until all cities have been visited
					Time: O(n) Space: O(n)'''
			while unvisited_cities and time.time()-start_time < time_allowance:
				'''Finds the index in the unvisited_cities list(a list of indexes) of the least cost city to travel to.
						Time: O(n) Space: O(n)'''
				idx = min(range(len(unvisited_cities)), key=lambda i: from_city.costTo(cities[unvisited_cities[i]]))
				'''Gets the actual city to travel to, rather than it's index'''
				to_city = cities[unvisited_cities[idx]]
				'''Terminate the visit to the city if it is impossible (cost is infinity)'''
				if from_city.costTo(to_city) == np.inf: break
				'''Add the city to the route and remove it from the unvisited_cities'''
				route.append(unvisited_cities.pop(idx))
				'''Update the city to start the next search from'''
				from_city = to_city
			'''We have found a possible solution'''
			count += 1
			'''If we have visitied all the cities, and there is a path back, update and return the found solution'''
			if not unvisited_cities and cities[route[-1]].costTo(cities[route[0]]) != np.inf:
				'''Create a TSPSolution form the converted route(list of indexes) to actual cities'''
				bssf = TSPSolution([cities[i] for i in route])
				break
		end_time = time.time()
		'''Return values of results'''
		results['cost'] = bssf.cost if bssf is not None else np.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''

	class State:
		"""
		A class to hold each sub-state of the Branch and Bound Algorithm.
		Store a cost_matrix of costs, a route, a list of remaining cities, and the current cost.
		The cities themselves are not stored within the State class, just their indexes;
		this was to hopefully make it faster.

		Complexity:
			Init Time: O(n^2)
			Copy Time: O(n^2)
			Update Time: O(n^2)
			Space: O(n^2)
		"""
		def __init__(self, init_cities):
			"""
			Creates the initial State of the Branch and Bound Algorithm and reduces it to find the base cost.
			Complexity:
				Time: O(n^2)
				Space: O(n^2)
			:param init_cities: List of cities to form cost_matrix from
			"""
			self.num_cities = len(init_cities)
			''' Create a new cost_matrix, each cell is filled with cost to move from city[row] -> city[column]
				The diagonals are marked as infinity so cities do not return to themselves
					Time: O(n^2) Space: O(n^2) '''
			self.cost_matrix = \
				np.array([[np.inf if i == j else init_cities[i].costTo(init_cities[j]) \
						   for j in range(self.num_cities)] for i in range(self.num_cities)])
			''' Initialize other default values '''
			self.last_city_idx = 0
			self.route = [0]
			''' Create list of unvisited cites indexes with all but 0 '''
			self.unvisited_cities = [i for i in range(1, self.num_cities)]
			''' Reduce the cost matrix to zero out minimum values '''
			self.cost = self.reduce_cost_matrix()

		def update_route(self, idx):
			"""
			This function is used to update the route by updating the route of the State.
			This involves re-reducing the matrix, and adding it to it's cost
			Complexity:
				Time: O(n^2)
				Space: O(n)

			:param idx: Index of city to move to; matches the corresponding row/col index of the cost matrix
			"""
			''' Add the cost to travel to the city to '''
			self.cost += self.cost_matrix[self.last_city_idx][idx]
			''' Remove the return edge from the matrix '''
			self.cost_matrix[idx][self.last_city_idx] = np.inf
			''' Remove the column that was visited '''
			self.cost_matrix[:, idx] = np.inf
			''' Remove the row that was visited '''
			self.cost_matrix[self.last_city_idx] = np.inf
			''' update the last visited city '''
			self.last_city_idx = idx
			''' Add the idx to the route '''
			self.route.append(idx)
			self.unvisited_cities.remove(idx)
			''' Reduce the cost matrix back to zeroed values
			 		Time: O(n^2) Space: O(n^2) '''
			self.cost += self.reduce_cost_matrix()

		def reduce_cost_matrix(self):
			"""
			Reduces all the rows in the cost matrix by the minimum value, then reduces the columns.
			This reduction skips rows that have a 0 or infinity as their minimum
			Complexity:
				Time: O(n^2)
				Space: O(n)

			:return: The base cost used to travel to any available city in the matrx
			"""
			def reduction(axis):
				"""
				Small helper function, reduces every row/col in the matrix by the minimum value of the row.
				Complexity:
					Time: O(n^2)
					Space: O(n)

				:param axis: The axis to reduce on, according to numpy's axis system
				:return: The sum cost removed by the reduction
				"""
				''' Use numpy to find minimum value in all rows/cols 
						Time: O(n) Space: O(n)'''
				min_vec = np.amin(self.cost_matrix, axis=axis)
				reduce_cost = 0
				''' Go through every value in the minimum values vector
						Time: O(n) Space: O(1)'''
				for i, c in enumerate(min_vec):
					if c == 0 or c == np.inf: continue
					reduce_cost += c
					''' subtract the minimum value from the row/col '''
					if axis == 1: self.cost_matrix[i] -= c
					else: self.cost_matrix[:,i] -= c
				''' return the sum of all minimum values '''
				return reduce_cost
			''' Return cost of the row reduction, followed by col reduction '''
			return reduction(axis=1) + reduction(axis=0)

		def get_next_state(self, idx):
			"""
			Creates a copy of the state, then moves it to the next city.
			None is returned if the move is invalid
			Complexity:
				Time: O(n^2) # None case is O(1)
				Space: O(n^2) # None case is O(1)

			:param idx: Index of the city to move the copied state to
			:return: None if no path is possible, otherwise a new State with an updated cost_matrix and route
			"""
			''' Terminate and return None if the city is unreachable '''
			if self.cost_matrix[self.last_city_idx][idx] == np.inf: return
			''' Create a deepcopy of the state to modify '''
			next_state = copy.deepcopy(self)
			''' Update the new state to the next city and return it '''
			next_state.update_route(idx)
			return next_state

		def __lt__(self, other):  # heapq uses < operator to push and pop right one
			"""
			Built in less than comparison; used to sort the heapq heap.
			Complexity:
				Time: O(1)
				Space: O(1)
			:param other: The other state to compare against
			:return: Boolean value of whether the state is less than the other state
			"""
			''' Check if the number of cities left is identical, if it is compare the costs
				ones with fewer cities to visit are prioritized '''
			if len(self.unvisited_cities) == len(other.unvisited_cities): return self.cost < other.cost
			else: return len(self.unvisited_cities) < len(other.unvisited_cities)

	def branchAndBound( self, time_allowance=60.0 ):
		"""
		This solves the TSP problem by using a greedy approach with the Nearest Neighbor (NN) algorithm.
		It chooses a random city, and then finds the next city it has not visited by getting the cheapest cost.
		The normal approach to this only runs on 1 index, but mine searches up to all the nodes to find a solution
		Complexity:
			Time: O(n!n^2)
			Space: O(n!n^2)

		:param time_allowance: Maximum time allowed to generate a solution for.
		:return: A TSP path it one was found otherwise np.inf
		"""
		''' Built in Heap library '''
		import heapq as hq
		'''Initialize the cites from the scenario'''
		results, cities = {}, list(self._scenario.getCities())
		'''Initialize some of the data we need for results'''
		ncities, count, idx_route = len(cities), 0, False
		''' Initialize the first State of the Branch and Bound '''
		states = [TSPSolver.State(cities)]
		'''Initialize some of the data we need for results'''
		states_pruned, total_states, max_queue_size = 0, 0, 0
		'''Initialize the first best solution so far [bssf] using the greedy algorithm
				Time: O(n^3) Space: O(n) '''
		bssf = self.greedy(time_allowance)['soln']
		start_time = time.time()
		''' Iterate through until we have searched all states or until we run out of time
				Iterations: n! Time: O(n!n^2) Space(n!n^2) '''
		while states and time.time()-start_time < time_allowance:
			''' Pop off the first item in the heap '''
			curr_state = hq.heappop(states)
			''' Confirm the cost of the state is still better than the bssf '''
			if curr_state.cost >= bssf.cost:
				states_pruned += 1
				continue
			''' Generate additional states for each non-visited cities 
				This Loop cycles through n-i states, one less for each time a state has been visited
					Time: O((n-i) n^2) Space: O((n-i) n^2) '''
			for idx in curr_state.unvisited_cities:
				''' Create a new State based off the current one with one additional city visited 
						Time: O(n^2) Space: O(n^2) '''
				new_state = curr_state.get_next_state(idx)
				total_states += 1
				if new_state is None: continue
				''' If the state has visited all the cities, see if it is a better solution than our bssf '''
				if len(new_state.route) == ncities:
					count += 1
					''' Update our bssf route if needed '''
					if new_state.cost < bssf.cost:
						bssf.route = new_state.route
						bssf.cost = new_state.cost
						idx_route = True
				elif new_state.cost < bssf.cost:
					''' Add the new states if they have a lower cost than the bssf '''
					hq.heappush(states, new_state)
				else:
					''' prune unnecessary states '''
					states_pruned += 1
			''' update max_queue_size for the results '''
			max_queue_size = max(max_queue_size, len(states))
		''' if we found a new route, we have to convert the list of indexes to a list of cities '''
		if idx_route: bssf.route = [cities[idx] for idx in bssf.route]
		end_time = time.time()
		'''Return values of results'''
		results['cost'] = bssf.cost if bssf is not None else np.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = max_queue_size
		results['total'] = total_states
		results['pruned'] = states_pruned
		return results



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		results, cities = {}, self._scenario.getCities()
		ncities, count, bssf = len(cities), 0, None
		start_time = time.time()
		bssf = self.greedy(time_allowance)['soln']
		if bssf.cost == np.inf: bssf = self.defaultRandomTour(time_allowance)['soln']
		cost = bssf.cost
		bssf.route = [c._index for c in bssf.route]
		og_cost = cost
		swap_count = 0
		def swap_cost(route, idxA, idxB):
			def getCosts(i):
				j = (i+1)%ncities
				return cities[route[i]].costTo(cities[route[j]]), cities[route[j]].costTo(cities[route[i]])
			012345 -> 043215
			cost: 0->1, 1->2, 2->3, 3->4, 4->5
			cost: 0->4, 4->3, 3->2, 2->1, 1->5
			if 0->4 and 1->5 == np.inf
				continue
			forward_cost = 0->4 + 4->5
			backward_cost = 0->4 + 1->5
			for off range(k+1): #k is 3
				idxA = (n+off)%ncities
				idxB = (n+off+1)%ncities
				forward_cost += idxA->idxB
				reverse_cost += indB->idxA
			if reverse_cost < forward_cost:
				preform_swap():
					A = firstCity
					while A->B:
						cityA.swap()
						if A == B: break
						A = A->prev
					end1 = cityA>next
					end2 = cityB>prev
					end1.next,A.next  = cityB, end2
					end2.prev,B.prev = cityA, end1
					cost += forward_cost-revere_cost

			reversed_cost = cities[route[idxA-1]].costTo(cities[route[idxB%ncities]]) + cities[route[(idxB+1)%ncities]].costTo(cities[route[idxA]])
			if reversed_cost == np.inf:
				return np.inf
			partial_cost = cities[route[idxA-1]].costTo(cities[route[idxA]]) + cities[route[idxB%ncities]].costTo(cities[route[(idxB+1)%ncities]])
			for i in range(idxA,idxB + 1):
				p,r = getCosts(i%ncities)
				partial_cost += p
				reversed_cost += r
				if reversed_cost == np.inf: return np.inf
			return reversed_cost - partial_cost

		def swap(route, cost, idxA, idxB):
			new_cost = swap_cost(route, idxA, idxB)
			if new_cost < 0:
				return True, route[:idxA] + route[idxB:idxA-1:-1] + route[idxB + 1%ncities:], cost+new_cost
			else:
				return False, route, cost
		k = 10
		is_updated = True
		while is_updated:
			is_updated = False
			if time.time()-start_time >= time_allowance: break
			for idx in range(ncities):
				if time.time()-start_time >= time_allowance: break
				for b in range(1, k + 1):
					change, route, cost = swap(route, cost, idx, idx+b)
					if change: swap_count += 1
					is_updated = is_updated or change
		bssf = TSPSolution([cities[idx] for idx in route])
		end_time = time.time()
		'''Return values of results'''
		results['cost'] = bssf.cost if bssf is not None else np.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = og_cost
		results['total'] = swap_count
		results['pruned'] = None
		return results



