import sys
import networkx as nx
import math
import matplotlib.pyplot as plt
from   scipy.spatial import Delaunay
import numpy as np

def euclidean(s, d):
    if( (s < 0) or (d < 0) or (s > len(P)) or (d > len(P)) ) :
        print("invalid query")
        sys.exit()
    xdiff = P[s][0]-P[d][0]
    ydiff = P[s][1]-P[d][1]
    return round(math.sqrt(xdiff**2+ydiff**2),3)
# ============== end of euclidean() ==============

G = nx.Graph()  # Graph
N = 0			# The number of point
w0 = 0			# Threshold value ( vertex to cycle )
P = []          # Point [ [x1, y1], [x2, y2], ...]  list
C = []          # Cycle vertex [ v1, v2, .. ] list
C_sz = 0        # The number of Cycle vertex
C_length = 0    # The Length of Cycle
C_valid = 1     # Cycle validation
C_max_dist = 0  # The longest distance (v)->(cycle)
C_max_v = 0     # longest distance vertex
C_max_E = []    # longest distance path
CE_IN = []      # Cycle edge in DG
CE_OUT = []     # Cycle edge out of DG

city_f = open("city-150.txt",'r')
cycle_f = open("test1.txt", 'r')


# Read cycle data and check
C_sz = int(cycle_f.readline())
C = list(map(int, cycle_f.readline().split()))

print(C[0])
print(C[-1])
print(len(C[1:]))
print(len(set(C[1:])))

if( C[0] != C[-1] or len(C[1:]) != len(set(C[1:])) ) :
    print("invalid cycle input")
    sys.exit()


# Read point data and check
N, w0 = map(int, city_f.readline().split())
for line in city_f :
    p = list(map(int, line.split()))
    if p in P :
        print( "invalid point input")
        sys.exit()
    else:
        P.append(p)


# Add node & e
# dge to G ( based on Delaunay Graph )
tri = Delaunay(P)
for i in range(1, N+1) :
	G.add_node( i, pos=P[i-1])
for t in tri.simplices.copy() :
    G.add_edge( t[0]+1, t[1]+1, weight=euclidean(t[0], t[1]))
    G.add_edge( t[1]+1, t[2]+1, weight=euclidean(t[1], t[2]))
    G.add_edge( t[0]+1, t[2]+1, weight=euclidean(t[0], t[2]))

# Cycle check
for i in range( len(C)-1 ) :
    e = [C[i], C[i+1]]
    if( e in G.edges() ): # cycle edge on DG
        CE_IN.append( e )
    else :
        CE_OUT.append( e )
    C_length += euclidean( C[i]-1, C[i+1]-1 )
if( len(CE_OUT) > 0 ) : C_valid = 0

# Calculate ditance between vertices and cycle
path = nx.johnson(G, weight='weight')   # or nx.all_pairs_dijkstra_path_length(G)
print("1. All Shortest Path (v-->cycle)")
print("vertex -- cycle_vertex ==> distance [path]")
for i in range( 1, len(P)+1 ) :
    if( i in C ) : continue
    else:
        min_dist = 1000*1000
        min_idx = 0
        for j in C :
            cur_dist = 0
            for k in range(len(path[i][j])-1):
                cur_dist += G[path[i][j][k]][path[i][j][k+1]]['weight']
            if cur_dist < min_dist:
                min_dist = cur_dist
                min_idx = j
        print( "   ", i, "--", min_idx, "==>", round(min_dist, 3), path[i][min_idx] )
        if C_max_dist < min_dist:
            C_max_dist = min_dist
            C_max_v = (i, min_idx)

if C_max_dist > w0:
    C_valid = 0
max_path = path[C_max_v[0]][C_max_v[1]]
for i in range(len(max_path)-1):
    C_max_E.append( [max_path[i], max_path[i+1]] )

# Report
print("\n 2. Longest Distance is", C_max_v, "==>", C_max_dist)
print("\n 3. Cycle Length is", C_length )
if( C_valid == 0 ) : print( "\n --> Wrong cycle answer!")
else :               print( "\n --> Correct cycle answer!")

G_subgraph = G.subgraph([29,7,30,16,58])
print(G_subgraph.edges())

# Draw Delaunay Graph and Cycle ========================================================================
pos = nx.get_node_attributes(G,'pos')
nx.draw_networkx_nodes(G, pos, node_size=100, alpha=1, node_color='skyblue')
# nx.draw_networkx_nodes (G, pos, node_size=300, nodelist=C, node_color='r')
# nx.draw_networkx_nodes (G, pos, node_size=300, nodelist=G_subgraph, node_color='r')
# nx.draw_networkx_edges (G, pos, edgelist=G_subgraph.edges(), width=3.0, alpha=0.6, style='dashed', edge_color='r') # Cycle edge on DG

nx.draw_networkx_nodes(G, pos, node_size=100, nodelist=[C_max_v[0]], node_color='orange')
nx.draw_networkx_edges(G, pos, width=0.5)
nx.draw_networkx_edges(G, pos, edgelist=CE_IN, width=7.0, alpha=1, edge_color='r') # Cycle edge on DG
nx.draw_networkx_edges(G, pos, edgelist=CE_OUT, width=3.0, alpha=0.6, style='dashed', edge_color='grey') # Cycle edge out of DG
nx.draw_networkx_edges(G, pos, edgelist=C_max_E, width=5.0, alpha=0.6, edge_color='orange') # longest distance
nx.draw_networkx_labels(G, pos, font_size=5, font_weight='bold')

edge_labels = nx.get_edge_attributes(G, 'weight')  # label로 표시한 변수를 선택
nx.draw_networkx_edge_labels(G, pos, font_size=3, edge_labels=edge_labels)

plt.tight_layout()
plt.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, hspace = 0, wspace = 0)
plt.show()
