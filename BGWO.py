import numpy as np
import random
import Connectivity_repair_SP as sp
import PropagationModel as pm
from copy import deepcopy
import networkx as nx
import time
import localSearch as ls

def metric_closure (G, weight='weight'):

    M = nx.Graph()
    Gnodes = set(G)

    # check for connected graph while processing first node
    all_paths_iter = nx.all_pairs_dijkstra(G, weight=weight)
    u, (distance, path) = next(all_paths_iter)
    if Gnodes - set(distance):
        msg = "G is not a connected graph. metric_closure is not defined."
        raise nx.NetworkXError(msg)
    Gnodes.remove(u)
    for v in Gnodes:
        M.add_edge(u, v, distance=distance[v], path=path[v])

    # first node done -- now process the rest
    for u, (distance, path) in all_paths_iter:
        Gnodes.remove(u)
        for v in Gnodes:
            M.add_edge(u, v, distance=distance[v], path=path[v])

    return M


class Obstacle:
    "type is a boolean variable. type=True if it is outer wall, False otherwise "
    def __init__(self, x0, y0, x1, y1,type):
        self.x0=x0
        self.x1 = x1
        self.y0=y0
        self.y1 = y1
        self.type=type


def load_obstacles(fichier):
    obstacles=[]
    with open(fichier, 'r') as file:
        for line in file:
            words = line.split(',')
            o=Obstacle(int(words[0]),int(words[1]),int(words[2]),int(words[3]),int(words[4])==0)
            obstacles.append(o)
    return obstacles


def generate_coordinates_archi1():

    coordinates=[]
    for j in range(5):
        for i in range(15):
            coordinates.append([i + 0.5, j+0.5])

    for j in range(2):
        for i in range(5):
            coordinates.append([i + 0.5, j+5.5])

    return coordinates



Rc = 10
Rs= 6
Ru=4
#nb_zones=20*35+40*30
nb_zones=15*5+5*2
pop_size= 25
epoch= 20
nb_targets =29
sensitivity= -96
file1 = open('targets_archi1', 'r')
list_target_points=[int(x) for x in file1.readline().split(',')]
print("time to create graph", end=" ")
obstacles= load_obstacles("archi1")
print("finish loading obstacles")
coordinates= generate_coordinates_archi1()
print("finish generating coordinates")
print("generating graph of zones....")
t= time.time()
graph = sp.generate_list_connections_between_positions(coordinates, obstacles,sensitivity, Rc)
print("finish generating graph of zones in ", end="  ")
print(time.time()-t)
print("generating metric closure....")
t= time.time()
metric = metric_closure(graph, weight='weight')
print("finish generating metric closure in ", end="  ")
print(time.time()-t)

def covering_zones_for_each_target():
    result=[]
    for target in list_target_points:
        lst=[]
        for i in range(len(coordinates)):
            if pm.Elfes_model(coordinates[i][0],coordinates[i][1],coordinates[target][0], coordinates[target][1], Rs, Ru):
                lst.append(i)
        result.append(lst)
    return result

target_covering_zones= covering_zones_for_each_target()


def create_random_pos():
    solution = [0, ] * nb_zones
    for i in range(nb_zones):
        if random.uniform(0, 1) > 0.75:
            solution[i] = 1
    return np.array(solution)

def calculate_cost(solution):
    return np.count_nonzero(solution)

def covered(target_point, solution, after_amend):
    #we do not need  the whole list of sensors covereing a target before the amend function bcause we will not compute augmented graph
    sensors_covering_targt = []
    for i in range(len(solution)):
        if solution[i] == 1 and i in target_covering_zones[target_point]:
            sensors_covering_targt.append(i)
            if not after_amend:
                return sensors_covering_targt

    return sensors_covering_targt

def calculate_coverage(wolf, sub_metric, augmented_graph, after_amend):
    #after_amend is added to cotrole either or not we augment agmented graph, this graph we use it only for
    #local search so there is no need to compute it before that.
    coverage = 0
    wolf.covered_targets=[]
    for i in range(len(list_target_points)):
        l = covered(i, wolf.position, after_amend)
        if len(l)>0:
            coverage = coverage + 1
            wolf.covered_targets.append(list_target_points[i])
            if after_amend:
                for k in l:
                    if not augmented_graph.has_edge(list_target_points[i],k):
                        augmented_graph.add_edge(list_target_points[i], k, weight=0)
                        sub_metric.add_edge(list_target_points[i], k, weight=0)

    wolf.coverage=coverage
    return coverage, augmented_graph, sub_metric


def compute_fitness(w, sub_metric, augmented_graph, after_amend, for_local_search):
    # we recompute coverage because after adding new sensors in amend function, other targets could be covered also. Yet we do not
    # recompute it after local search because we know that the coverage will to be neither decreased nor increased.
    if for_local_search:
        coverage = w.coverage
    else:
        coverage,  augmented_graph, sub_metric = calculate_coverage(w, sub_metric, augmented_graph, after_amend)
    cost = calculate_cost(w.position)
    w.fitness = (nb_targets - coverage) / nb_targets + 0.5*(cost / nb_zones)
    return augmented_graph, sub_metric


class wolf:
    def __init__(self,):
        self.position = []
        self.fitness = []
        self.coverage=0
        self.covered_targets=[]


def sig(x):
 return 1/(1 + np.exp(-10*x-0.5))


def sigmoid_transformation(position):
    pos=[]
    for i in range(nb_zones):
        r=random.random()
        s= sig(position[i])
        if s < r:
            pos.append(1)
        else:
            pos.append(0)
    return np.array(pos)


def get_best_solutions(pop, best=3):
    pop = sorted(pop, key=lambda agent: agent.fitness)
    return deepcopy(pop[:best])


def create_pop():
    pop=[]
    for i in range(pop_size):
        w= wolf()
        w.position= create_random_pos()
        g = compute_fitness(w,None,None, False, False)
        pop.append(w)
    return pop


def create_deployment_graph(solution):
    index_of_deployed_sensors = [i for i in range(len(solution)) if solution[i] == 1]
    deployment_graph = graph.subgraph(index_of_deployed_sensors).copy()
    return deployment_graph


def amend_position(solution):
    solution = sigmoid_transformation(solution)
    deployment_graph= create_deployment_graph(solution)
    disjoint_sets = sp.distinct_connected_components(deployment_graph)
    if len(disjoint_sets) > 1:
        sp.connectivity_repair_heuristic(disjoint_sets, solution,  graph, coordinates)
    return np.array(solution), deployment_graph


def solve():

    pop= create_pop()
    for current_epoch in range(0, epoch):

        a = 2 - 2 * current_epoch / (epoch - 1)
        list_best = get_best_solutions(pop, best=3)
        pop_new = []
        for idx in range(0, pop_size):
            A1, A2, A3 = a * (2 * np.random.uniform() - 1), a * (2 * np.random.uniform() - 1), a * (2 * np.random.uniform() - 1)
            C1, C2, C3 = 2 * np.random.uniform(), 2 * np.random.uniform(), 2 * np.random.uniform()
            X1 = np.abs(list_best[0].position - A1 * np.abs(C1 * list_best[0].position - pop[idx].position))
            X2 = np.abs(list_best[1].position - A2 * np.abs(C2 * list_best[1].position - pop[idx].position))
            X3 = np.abs(list_best[2].position - A3 * np.abs(C3 * list_best[2].position - pop[idx].position))
            pos_new = (X1 + X2 + X3) / 3.0
            pos_new, graph = amend_position(pos_new)
            w = wolf()
            w.position = pos_new
            sub_metric=metric.subgraph(list(graph.nodes)).copy()
            augmented_graph = graph.subgraph(list(graph.nodes)).copy()
            """
            print("augmented graph before fitness  ",end="  ")
            print(list(augmented_graph.edges()))
            """
            augmented_graph,sub_metric = compute_fitness(w, sub_metric, augmented_graph, True, False)
            """
            print("augmented graph after fitness  ", end="  ")
            print(list(augmented_graph.edges()))
            print("list of covered point before local search ", end="  ")
            print(w.covered_targets)
            """
            """
            print("*******************************************")
            print("nb nodes in graph ", end=" ")
            print(len(graph.nodes), end="  ")
            print("nb target covered ")
            print(len(w.covered_targets), end="  ")
            print(" nb nodes in augmented nodes ", end="  ")
            print(len(augmented_graph.nodes))
            print("*******************************************")
            """
            print("***********************************************************************")
            coverage_before_ls, g,c = calculate_coverage(w,None, None, False)
            n=[]
            for i in range(len(w.position)):
                if w.position[i]==1:
                    n.append(i)
            print("check if the sensors really cover all the target points")
            print("list of targets : ", end="  ")
            print(list_target_points)
            print("deployed sensors", end=" ")
            print(n)
            remove_nodes = ls.local_search(sub_metric, graph, augmented_graph, metric, w.covered_targets)
            for r in remove_nodes:
                w.position[r]=0
            n=[]
            for i in range(len(w.position)):
                if w.position[i]==1:
                    n.append(i)
            print("deployed sensors after removing", end=" ")
            print(n)
            g = compute_fitness(w, augmented_graph, False, True)
            graph = create_deployment_graph(w.position)
            print("after local search : ", end="  ")
            print("connected? ", end=" ")
            s, g = calculate_coverage(w, None, False)
            print(len(sp.distinct_connected_components(graph)) == 1, end=" ")
            print("modified coverage ? ", end="  ")
            #print(w.coverage, end=" ")
            print(s)
            print(coverage_before_ls)
            #print(w.coverage == ms)
           #print(s == ms)
            pop_new.append(deepcopy(w))
        pop_new.extend(deepcopy(list_best))
        pop_new = sorted(pop_new, key=lambda agent: agent.fitness)
        pop=deepcopy(pop_new[:3])
        pop.extend(deepcopy(pop_new[3:pop_size]))
    best=get_best_solutions(pop, best=1)
    print("best sol ", end="  ")
    print(best[0].position)
    print(best[0].fitness)
    print("check if the sensors really cover all the target points")
    print("list of targets : ", end="  ")
    print(list_target_points)
    print("coverage of best sol ", end="  ")
    print(best[0].coverage)
    for i in range(len(best[0].position)):
        if best[0].position[i]==1:
            print("sensor at position ", end=" ")
            print(i)
            for j in range(len(list_target_points)):
                if i in target_covering_zones[j]:
                    print("sensor ", end=" ")
                    print(i, end=" ")
                    print("covers ", end=" ")
                    print(j)

    return best
solve()