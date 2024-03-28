import networkx as nx
import networkx.algorithms.components as comp
import math
import PropagationModel as pm

# input must be a graph : either the deployment graph or the tree in SP
def distinct_connected_components(graph):
    connected_nodes = []
    for connected_set in list(comp.connected_components(graph)):
        connected_nodes.append(list(connected_set))

    return connected_nodes


# node1 and node2 are indexes, list_disjoints_neighbors is list of list
def in_same_set(node1, node2, list_of_list):
    for l in list(list_of_list):
        if node1 in list(l):
            if node2 in list(l):
                return True
            else:
                return False
def distance(ax, ay, bx, by):
    return math.sqrt(pow(ax - bx, 2) + pow(ay - by, 2))


def second_step_of_connectivity_repair(list_disjoints_neighbors, graph, zones):
    terminals = []
    tree = nx.Graph()
    for s in list_disjoints_neighbors:
        terminals.append(s[0])
    num_terminals = len(terminals)
    distance_between_terminals = []
    for i in range(num_terminals):
        for j in range(i + 1, num_terminals):
            pair = dict()
            pair['src'] = terminals[i]
            pair['dest'] = terminals[j]
            pair['dis'] = distance(zones[terminals[i]].zonepolygon.centroid.x, zones[terminals[i]].zonepolygon.centroid.y,
                                          zones[terminals[j]].zonepolygon.centroid.x, zones[terminals[j]].zonepolygon.centroid.y)
            pair['inserted'] = False
            distance_between_terminals.append(pair)
    distance_between_terminals = sorted(distance_between_terminals, key=lambda d: d['dis'])
    inserted_terminals = []
    i = 0
    while len(inserted_terminals) < len(terminals):
        if distance_between_terminals[i]['src'] not in inserted_terminals or distance_between_terminals[i]['dest'] not in inserted_terminals:
            least_cost_path = nx.shortest_path(graph, source=distance_between_terminals[i]['src'],
                                               target=distance_between_terminals[i]['dest'])
            nx.add_path(tree, least_cost_path)
            inserted_terminals.append(distance_between_terminals[i]['src'])
            inserted_terminals.append(distance_between_terminals[i]['dest'])
            inserted_terminals = list(set(inserted_terminals))
            distance_between_terminals[i]['inserted'] = True

            # if i==len(distance_between_terminals):
            #    break
        i = i + 1

    D = distinct_connected_components(tree)
    if len(D) > 1:
        not_added_paths = [elem for elem in distance_between_terminals if elem['inserted'] == False]
        f = 0
        while len(D) > 1:
            elem = not_added_paths[f]
            f = f + 1
            if in_same_set(elem['src'], elem['dest'], D) == False:
                least_cost_path = nx.shortest_path(graph, source=elem['src'], target=elem['dest'])
                nx.add_path(tree, least_cost_path)
                inserted_terminals.append(elem['src'])
                inserted_terminals.append(elem['dest'])
                D = distinct_connected_components(tree)

    removed_nodes = []
    """
    borders = [node for node in list(tree.nodes) if tree.degree(node) == 1]
  
    for b in borders:
        e = tree.edges(b)
        e = list(list(e)[0])
        removed_node = b
        if removed_node == e[0]:
            k = e[1]
        else:
            k = [0]

        while in_same_set(e[0], e[1], list_disjoints_neighbors):
            # removed_nodes.append(b)
            tree.remove_node(removed_node)
            if tree.degree(k) == 1:
                e = tree.edges(k)
                e = list(list(e)[0])
                removed_node = k
                if removed_node == e[0]:
                    k = e[1]
                else:
                    k = [0]

            else:
                break
    """
    return tree, removed_nodes


def get_neighbors_positions_of_disjoints_connected_set(lst_disjoints_connected_set, lst_neighbors_per_zone):
    S = [[] for x in range(len(lst_disjoints_connected_set))]
    i = 0
    for d in lst_disjoints_connected_set:
        s = []
        for p in d:
            s.extend(lst_neighbors_per_zone[p])
            s = [elem for elem in s if elem not in d]
        S[i] = list(set(s))
        i = i + 1
    return S


def highest_occurence(lst):
    return max(lst, key=lst.count)


def index_of(a, lst):
    for i in range(0, len(lst)):
        if lst[i] == a:
            return i
    return -1


def first_step_of_connectivity_repair(individual, disjoint_sets, lst_neighbors_per_zone, lst_neighbors_per_set):
    cont = True
    while cont:
        P = []
        for s in lst_neighbors_per_set:
            P = P + s
        l = highest_occurence(P)
        if l > 1:
            individual.deployment[l] = 1
            merged_set = []
            merged_neighbors = []
            lis_index = []
            for i in range(0, len(lst_neighbors_per_set)):
                ind = index_of(l, lst_neighbors_per_set[i])
                if ind != -1:
                    lis_index.append(i)
            for k in lis_index:
                merged_set = merged_set + disjoint_sets[k]
                merged_neighbors.extend((lst_neighbors_per_set[k]))

            for k in reversed(lis_index):
                disjoint_sets.pop(k)
                lst_neighbors_per_set.pop(k)

            merged_set.append(l)
            disjoint_sets.append(merged_set)
            merged_neighbors.extend(lst_neighbors_per_zone[l])
            merged_neighbors = list(set(merged_neighbors))
            merged_neighbors.remove(l)
            lst_neighbors_per_set.append(merged_neighbors)
        else:
            cont = False


def generate_list_connections_between_positions(zones, obstacles, sensitivity, Rc):
    graph = nx.Graph()
    for i in range(len(zones)-1):
        for j in range(i + 1, len(zones)):
            if check_connectivity(zones[i].zonepolygon.centroid.x, zones[i].zonepolygon.centroid.y, zones[j].zonepolygon.centroid.x, zones[j].zonepolygon.centroid.y, obstacles,
                                  sensitivity, Rc):
                if not graph.has_edge(i, j):
                    graph.add_edge(i, j)

    return graph

def generate_coverage_graph(zones, Rs, Ru):
    graph = nx.Graph()
    for i in range(len(zones)-1):
        for j in range(i + 1, len(zones)):
            if pm.Elfes_model(zones[i].zonepolygon.centroid.x, zones[i].zonepolygon.centroid.y, zones[j].zonepolygon.centroid.x, zones[j].zonepolygon.centroid.y, Rs,
                                  Ru):
                if not graph.has_edge(i, j):
                    graph.add_edge(i, j)

    return graph



def check_connectivity(position_1x, position_1y, position_2x, position_2y, obstacles, threshold, rc):
    return pm.MWM(position_1x, position_1y, position_2x, position_2y, obstacles, threshold, rc)

from copy import deepcopy

def connectivity_repair_heuristic(disjoint_sets, individual, graph, zones):
    steiner_tree, removed_nodes = second_step_of_connectivity_repair(disjoint_sets, graph, zones)
    for i in list(set(steiner_tree.nodes)):
        individual[i] = 1
    return deepcopy(individual)

