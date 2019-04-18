import sys
import networkx as nx
import math
import matplotlib.pyplot as plt
from   scipy.spatial import Delaunay
import numpy as np
import copy
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


def check_cycle_length(G, C, P):
    C_length = 0
    edge_in = []
    edge_out = []
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

    if C_max_dist > w0:
        return False
    else:
        return True


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


def make_cycle(G, V):
    avg_x = 0
    avg_y = 0
    for v in V:
        avg_x += G.node[v]['pos'][0]
        avg_y += G.node[v]['pos'][1]
    avg_x = avg_x / len(V)
    avg_y = avg_y / len(V)

    # Sort the vertexes by using radius
    radis = []
    for v in V:
        x = G.node[v]['pos'][0]- avg_x
        y = G.node[v]['pos'][1]- avg_y
        rad = np.arctan2(y, x)
        if rad < 0:
            rad += 2*math.pi
        radis.append(rad)

    # Make cycle from pins
    v_rads = dict(zip(V, radis))
    sorted_pins = sorted(v_rads.items(), key=lambda k : k[1])
    cycle = []
    G_temp = copy.deepcopy(G)    ## Copy the graph G
    for i in range(len(sorted_pins) - 1):
        line_path = nx.shortest_path(G_temp, sorted_pins[i][0], sorted_pins[i+1][0], weight='weight')
        cycle = cycle + line_path
        G_temp.remove_nodes_from(line_path[1:len(line_path)-1])
    line_path = nx.shortest_path(G_temp, sorted_pins[len(sorted_pins) - 1][0], sorted_pins[0][0], weight='weight')
    cycle = cycle + line_path
    return cycle


def target_function(G, P, w, boundary, **pins):
    V = select_vertexes_from_points(G, **pins)
    cycle = make_cycle(G, V)

    print(cycle)
    C_length , edge_in, edge_out = check_cycle_length(G, cycle, P)

    #graph_vis(G, cycle, V, edge_in, edge_out)

    worst_score = ((boundary[2]-boundary[0])+(boundary[3]-boundary[1]))*2

    if not constraint_correct(G, C, P, w):
        return -1*worst_score
    else:
        return -1*C_length


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


def graph_vis(G, cycle_v=None, pins=None, CE_IN=None, CE_OUT=None):
    pos = nx.get_node_attributes(G,'pos')

    nx.draw_networkx_nodes(G, pos, node_size=300, alpha=0.7, node_color='skyblue')
    if not cycle_v ==None:
        nx.draw_networkx_nodes(G, pos, node_size=300, nodelist=cycle_v, node_color='r')
        nx.draw_networkx_nodes(G, pos, node_size=300, nodelist=pins, node_color='y')
    nx.draw_networkx_edges (G, pos, width=0.5)
    if not CE_IN == None:
        nx.draw_networkx_edges(G, pos, edgelist=CE_IN, width=3.0, alpha=0.6, style='dashed', edge_color='r',
                                label=nx.get_edge_attributes(G, 'weight')) # Cycle edge out of DG
    if not CE_OUT == None:
        nx.draw_networkx_edges(G, pos, edgelist=CE_OUT, width=3.0, alpha=0.6, style='dashed', edge_color='grey',
                                label=nx.get_edge_attributes(G, 'weight')) # Cycle edge out of DG

    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold' )

    plt.tight_layout()
    plt.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, hspace = 0, wspace = 0)
    plt.show()


C = read_cycle("citycycle-060.txt")
P, threshold = read_points("city-060.txt")
G = makeDelaunay(P)
C_length, CE_IN, CE_OUT = check_cycle_length(G, C, P)
min_x, min_y, max_x, max_y = boundary(P)
min_val = min(min_x, min_y)
max_val = max(max_x, max_y)

print(" Threshold : {}".format(threshold))

# Make an argument of target function be only cycle vertexes
vertex_n =8
bound = [(min_val, max_val) for x in range(vertex_n)]
pbounds = dict(zip(map(str,range(vertex_n+1)), bound))
func = lambda **pins : target_function(G, P, threshold, (min_x, min_y, max_x, max_y), **pins)
optimizer = BayesianOptimization(
    f = func,
    pbounds=pbounds,
    verbose=2
)
optimizer.maximize(init_points=vertex_n, acq='ei', n_iter=100)

print(optimizer.max)
# Report
print( "\n 3. Cycle Length is", optimizer.max['target'] )

# Draw Delaunay Graph and Cycle ========================================================================
