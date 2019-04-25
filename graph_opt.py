import sys
import networkx as nx
import math
import matplotlib.pyplot as plt
from   scipy.spatial import Delaunay
import numpy as np
import copy
import random
from deap import base, creator, tools
import multiprocessing

from bayes_opt import BayesianOptimization

def euclidean( P, s, d ) :
    if( (s < 0) or (d < 0) or (s > len(P)) or (d > len(P)) ) :
        print("invalid query")
        sys.exit()
    xdiff = P[s][0]-P[d][0]
    ydiff = P[s][1]-P[d][1]
    return round(math.sqrt(xdiff**2+ydiff**2),3)


def read_cycle(file_name):
    cycle_f = open(file_name, 'r')
    cycle_f.readline()
    C = list(map(int, cycle_f.readline().split()))
    if( C[0] != C[-1] or len(C[1:]) != len(set(C[1:])) ) :
        print("invalid cycle input")
        sys.exit()
    return C


def read_points(file_name):
    city_f = open(file_name,'r')
    N, treshold = map(int, city_f.readline().split())
    P = []
    for line in city_f:
        p = list(map(int, line.split()))
        if( p in P ) :
            print( "invalid point input")
            sys.exit()
        else:
            P.append( p )
    return P, treshold


def write_cycle(cycle, file_name):
    cycle_f = open(file_name, 'w')
    cycle_f.write(str(len(cycle)) + '\n')
    for c in cycle:
        cycle_f.write(str(c) + " ")
    cycle_f.close()


def makeDelaunay(P):
    G = nx.Graph()
    tri = Delaunay(P)
    N = len(P)
    for i in range(1, N+1) :
        G.add_node( i, pos=P[i-1])
    for t in tri.simplices.copy() :
        G.add_edge( t[0]+1, t[1]+1, weight=euclidean(P, t[0], t[1]))
        G.add_edge( t[1]+1, t[2]+1, weight=euclidean(P, t[1], t[2]))
        G.add_edge( t[0]+1, t[2]+1, weight=euclidean(P, t[0], t[2]))
    return G


def remove_duplicate_in_cycle(C):
    before = C[0]
    cycle = [C[0]]
    for c_i in C:
        if c_i == before:
            continue
        else:
            before = c_i
            cycle.append(c_i)
    return cycle


def check_cycle_length(G, C, P):
    C_length = 0
    edge_in = []
    edge_out = []
    j = 0
    while j < len(C) - 1:
        if C[j] == C[j+1]:
            C.pop(j)
        j += 1

    for i in range(len(C)-1):
        e = [C[i], C[i+1]]
        if( e in G.edges() ): # cycle edge on DG
            edge_in.append( e )
        else :
            edge_out.append( e )
        C_length += euclidean(P, C[i]-1, C[i+1]-1 )

    return C_length, edge_in, edge_out


def constraint_correct(G, C, P, w0):
    path = nx.johnson(G, weight='weight')   # or nx.all_pairs_dijkstra_path_length(G)
    C_max_dist = 0
    for i in range( 1, len(P)+1 ) :
        if( i in C ) : continue
        else:
            min_dist = 1000*1000
            min_idx = 0
            for j in C :
                cur_dist = 0
                for k in range( len(path[i][j])-1 ) :
                    cur_dist += G[path[i][j][k]][path[i][j][k+1]]['weight']
                if( cur_dist < min_dist ) :
                    min_dist = cur_dist
                    min_idx = j
            if( C_max_dist < min_dist ):
                C_max_dist = min_dist
                C_max_v = (i, min_idx)

    max_path = path[C_max_v[0]][C_max_v[1]]
    max_edges = []
    for i in range(len(max_path)-1):
        max_edges.append([max_path[i], max_path[i+1]])

    if C_max_dist > w0:
        return False, max_edges
    else:
        return True, max_edges


def select_vertexes_from_points(G, **pins):
    V = [int(x) for x in pins.values()]

    i = 0
    points = []
    for j in range(int(len(V)/2)):
        points.append((V[i],V[i+1]))
        i += 2

    vertexes = []
    for p in points:
        closest_node = 1
        for n in G.nodes:
            distance = lambda x, y : math.sqrt((x['pos'][0] - y[0])**2 + (x['pos'][1] - y[1])**2)
            if distance(G.node[closest_node], p) > distance(G.node[n], p):
                closest_node = n
        vertexes.append(closest_node)
    return vertexes


def make_cycle(G, V, center, avg_center=False):
    if avg_center == True:
        avg_x = 0
        avg_y = 0
        for v in V:
            avg_x += G.node[v]['pos'][0]
            avg_y += G.node[v]['pos'][1]
        avg_x = avg_x / len(V)
        avg_y = avg_y / len(V)
        center = [avg_x, avg_y]
    # Sort the vertexes by using radius
    radis = []
    for v in V:
        x = G.node[v]['pos'][0]- center[0]
        y = G.node[v]['pos'][1]- center[1]
        rad = np.arctan2(y, x)
        if rad < 0:
            rad += 2*math.pi
        radis.append(rad)

    # Make cycle from pins
    v_rads = dict(zip(V, radis))
    sorted_pins = sorted(v_rads.items(), key=lambda k : k[1])
    pins = []
    for i in sorted_pins:
        pins.append(i[0])
    cycle = []
    G_temp = copy.deepcopy(G)    ## Copy the graph G
    i = 0
    while i < (len(pins) - 1):
        line_path = nx.shortest_path(G_temp, pins[i], pins[i+1], weight='weight')
        cycle = cycle + line_path
        intersect = list(set(line_path[1:len(line_path)-1]) & set(pins))
        if len(intersect) == 0:
            G_temp.remove_nodes_from(line_path[1:len(line_path)-1])
            i += 1
        else:
            G_temp.remove_nodes_from(line_path[1:len(line_path)-1])
            for inter_ in intersect:
                j = 0
                while j < len(pins):
                    if inter_ == pins[j]:
                        pins.pop(j)
                    else:
                        j += 1
            i += 1
    line_path = nx.shortest_path(G_temp, pins[len(pins) - 1], pins[0], weight='weight')
    cycle = cycle + line_path
    cycle = remove_duplicate_in_cycle(cycle)
    return cycle


def target_function(G, P, w, boundary, case=None, **pins):
    if case == 'geo':
        V = select_vertexes_from_points(G, **pins)
    elif case == 'int':
        V = [int(x) for x in pins.values()]
    else:
        exit()
    worst_score = ((boundary[2]-boundary[0])+(boundary[3]-boundary[1]))*2
    try:
        cycle = make_cycle(G, V)
    except:
        print("|Cycle build failed 1 |", end="")
        return -1
    try:
        C_length , edge_in, edge_out = check_cycle_length(G, cycle, P)
    except:
        print("|Cycle build failed 2 |", end="")
        return -1

    is_in, path = constraint_correct(G,cycle,P,w)
    #graph_vis(G, cycle, V, edge_in, edge_out, path)

    if not is_in:
        print("|Constraint failed    |", end="")
        return -1
    else:
        print("|Cycle :: {:.2f}     |".format(C_length), end="")
        return -C_length/worst_score


def ga_evaluation(G, P, w, center, V):
    worst_score = G.size(weight='weight')

    try:
        cycle = make_cycle(G, V, center)
    except:
        return [worst_score]
    try:
        C_length , edge_in, edge_out = check_cycle_length(G, cycle, P)
    except:
        return [worst_score-1]

    is_in, path = constraint_correct(G,cycle,P,w)
    #graph_vis(G, cycle, V, edge_in, edge_out, path)

    if not is_in:
        return [worst_score]
    else:
        return [C_length]


def boundary(P):
    min_x = 100000
    max_x = -100000
    min_y = 100000
    max_y = -100000
    for i in P:
        if min_x > i[0]:
            min_x = i[0]
        if max_x < i[0]:
            max_x = i[0]
        if min_y > i[1]:
            min_y = i[1]
        if max_y < i[1]:
            max_y = i[1]
    return min_x, min_y, max_x, max_y


def graph_vis(G, cycle_v=None, pins=None, CE_IN=None, CE_OUT=None, path=None):
    pos = nx.get_node_attributes(G,'pos')
    nx.draw_networkx_nodes(G, pos, node_size=300, alpha=0.7, node_color='skyblue')
    if not cycle_v ==None:
        nx.draw_networkx_nodes(G, pos, node_size=300, nodelist=cycle_v, node_color='r')
    if not pins == None:
        nx.draw_networkx_nodes(G, pos, node_size=300, nodelist=pins, node_color='y')

    nx.draw_networkx_edges (G, pos, width=0.5)
    if not CE_IN == None:
        nx.draw_networkx_edges(G, pos, edgelist=CE_IN, width=3.0, alpha=0.6, style='dashed', edge_color='r',
                                label=nx.get_edge_attributes(G, 'weight')) # Cycle edge out of DG
    if not CE_OUT == None:
        nx.draw_networkx_edges(G, pos, edgelist=CE_OUT, width=3.0, alpha=0.6, style='dashed', edge_color='grey',
                                label=nx.get_edge_attributes(G, 'weight')) # Cycle edge out of DG
    if not path == None:
        nx.draw_networkx_edges(G, pos, edgelist=path, width=3.0, alpha=0.6, style='dashed', edge_color='y',
                               label=nx.get_edge_attributes(G, 'weight')) # Cycle edge out of DG

    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold' )

    plt.tight_layout()
    plt.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, hspace = 0, wspace = 0)
    plt.show()


def genetic_optimization(toolbox, population_size=100, cross_prob=0.5, mutation_prob=0.1, generations=100):
    pop = toolbox.population(n=population_size)
    halloffame = tools.HallOfFame(maxsize=1)
    i = 0
    fitnesses = map(toolbox.evaluate, pop)
    for ind, fit in zip(pop, fitnesses):
        print("Pop {} Init".format(i))
        ind.fitness.values = fit
        i += 1

    for g in range(generations):
        print("Generation {}: ".format(g+1), end="")
        next_gen = toolbox.select(pop, int(population_size/10), int(population_size/3))
        next_gen = map(toolbox.clone, next_gen)
        next_gen = list(next_gen)
        #Crossover
        for child1, child2 in zip(next_gen[::2], next_gen[1::2]):
            if random.random() < cross_prob:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values

        #Mutation
        for mutant in next_gen:
            if random.random() < mutation_prob:
                toolbox.mutate(mutant)
                del mutant.fitness.values

        invalid_ind = [ind for ind in next_gen if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        pop[:] = next_gen
        fits = [ind.fitness.values[0] for ind in pop]
        halloffame.update(pop)
        print("AVG:{:.2f}  MIN:{:.2f}  MAX:{:.2f} {}".format(np.mean(fits), np.min(fits), np.max(fits), halloffame[0]))

    return halloffame[0]

### 0. Init values
P, threshold = read_points("city-200.txt")
G = makeDelaunay(P)
min_x, min_y, max_x, max_y = boundary(P)
min_val = min(min_x, min_y)
max_val = max(max_x, max_y)
center = ((min_x+max_x)/2, (min_y+max_y)/2)
print(" Threshold : {}, Size of graph : {}".format(threshold, len(G)))

### 1. Gaussian Process :: ( Bayesian Optimization ) > It takes too long because kernel function calculations.
### If you have to use GP, it can apply to adjust parameter of GA especailly "pin_num" which reperesent
### the number of pin set.
"""
vertex_n = 4
bound = [(min_val, max_val) for x in range(vertex_n)]
bound_integer = [(1, len(G)+1) for x in range(vertex_n)]
pbounds = dict(zip(map(str,range(vertex_n+1)), bound))
pbounds_int = dict(zip(map(str,range(vertex_n+1)), bound_integer))
func = lambda **pins : target_function(G, P, threshold, (min_x, min_y, max_x, max_y), case="geo", **pins)
func_int = lambda **pins : target_function(G, P, threshold, (min_x, min_y, max_x, max_y), case="int", **pins)
optimizer = BayesianOptimization(
    f = func_int,
    pbounds=pbounds_int,
    verbose=2
)
optimizer.maximize(init_points=100, acq='ei', n_iter=1000)
"""

### 2. Genetic Algorithm as an integer programming solver
pin_num = 30
pop_size = 300
#pool = multiprocessing.Pool()
toolbox = base.Toolbox()
#toolbox.register("map", pool.map)
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)
toolbox.register("indices", random.sample, range(pin_num), len(G))
toolbox.register("attribute", lambda :random.randint(1, len(G)))
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attribute, n=pin_num)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("mate", tools.cxUniform, indpb=2/pin_num)
toolbox.register("mutate", tools.mutUniformInt, low=1, up=len(P), indpb=3/pin_num)
toolbox.register("select", tools.selTournament)
toolbox.register("evaluate", lambda pins : ga_evaluation(G, P, threshold, center, pins))

B_E_S_T = genetic_optimization(toolbox, generations=500, population_size=pop_size, mutation_prob=0.5, cross_prob=0.3)
cycle = make_cycle(G, B_E_S_T, center)
constraint, far_edges = constraint_correct(G, cycle, P, threshold)
length, edge_in, edge_out = check_cycle_length(G, cycle, P)
graph_vis(G, cycle, B_E_S_T, edge_in, edge_out, far_edges)

write_cycle(cycle, "test1.txt")
