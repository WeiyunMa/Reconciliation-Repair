# python libraries
from cStringIO import StringIO

# BioPython libraries
from Bio import Phylo

import DP
import pickle
import networkx as nx
from ete2 import Tree
import newickFormatReader
from ReconciliationGraph import buildReconciliation
import MasterReconciliation
import os.path

# def find_root(tree):
#     nodes = tree.keys()
#     for node in tree.keys():
#         _, child1, child2 = tree[node]
#         if child1 in nodes:
#             nodes.remove(child1)
#         if child2 in nodes:
#             nodes.remove(child2)
#     return nodes[0]

# def edge_to_vertex_tree(tree):
#     nodes = {}
#     for edge in tree:
#         _, current_node, child_edge_1, child_edge_2 = tree[edge]
#         if child_edge_1 is not None:
#             _, child_node_1 = child_edge_1
#         else:
#             child_node_1 = None
#         if child_edge_2 is not None:
#             _, child_node_2 = child_edge_2
#         else:
#             child_node_2 = None
#         nodes[current_node] = (child_node_1, child_node_2)
#     return nodes


# def get_dist(S, s1, s2):
#    dist = 0
#    while S.search_nodes(name=s2)[0].name != s1 and S.get_common_ancestor(s1, s2).name != s2:
#        s2 = S.search_nodes(name=s2)[0].up.name
#        dist += 1
#     return dist
    
    
def recon_tree_to_dtl(T):
    sigma, delta, theta, xi = [], [], [], []
    M, L, tau = {}, {}, {}
    
    for mapping_node in T.keys():
        g, s = mapping_node
        # if g in M.keys():
        #     if S.get_common_ancestor(M[g], s).name == M[g]: 
        #         M[g] = s 
        # else:
        #     M[g] = s

        event, _, _ = T[mapping_node]
        if event != 'L':
            M[g] = s
            if event == 'S':
                sigma.append(g)
            elif event == 'D':
                delta.append(g)
            elif event == 'C':
                L[g] = s
            elif event == 'T':
                theta.append(g)
            
    for g in theta:
        _, child_1, child_2 = T[(g, M[g])]
        g1, s1 = child_1
        g2, s2 = child_2
        # if child_1 == (None, None) or child_2 == (None, None):
        #    print "what", (g, M[g])
        if M[g] == s1: # M[g1] or S.get_common_ancestor(M[g], M[g1]).name == M[g]:
            # print "yay"
            # print (g, g2)
            # print s2
            xi.append((g, g2))
            tau[g] = s2
        else:
            # print (g, g1)
            # print s1
            xi.append((g, g1))
            tau[g] = s1
            
    return (L, M, sigma, delta, theta, xi, tau)
        

def dtl_to_recon_tree(S, G, alpha):
    L, M, sigma, delta, theta, xi, tau = alpha
    T = {}
    # print S
    for g in [node.name for node in G.traverse("preorder")]:
        if g in sigma:
            g_value = ['S']
        elif g in theta:
            g_value = ['T']
        elif g in delta:
            g_value = ['D']
        else:
            g_value = ['C']
            
        for g_child in [node.name for node in G.search_nodes(name=g)[0].children]:
            if g in sigma:
                # print 'S', M[g], M[g_child]
                dist = 0 if M[g] == M[g_child] else int(S.get_distance(M[g], M[g_child]))
                pa = M[g]
            elif g in theta and (g, g_child) in xi:
                # print 'T', M[g], tau[g], M[g_child]
                dist = 1 if tau[g] == M[g_child] else int(S.get_distance(tau[g], M[g_child])) + 1
                pa = tau[g]
            else:
                # print 'D', M[g], M[g_child]
                dist = 1 if M[g] == M[g_child] else int(S.get_distance(M[g], M[g_child])) + 1
                pa = M[g]
            # print M[g], M[g_child]
            # print g, g_child
            # if g in sigma:
            #     # print 'S', M[g], M[g_child]
            #     dist = get_dist(S, M[g], M[g_child])
            # elif g in theta and (g, g_child) in xi:
            #     # print 'T', tau[g], M[g_child]
            #     dist = get_dist(S, tau[g], M[g_child]) + 1
            # else:
            #     # print 'D', M[g], M[g_child]
            #     dist = get_dist(S, M[g], M[g_child]) + 1
            
            x = M[g_child]
            if x == pa or S.search_nodes(name=x)[0] in S.search_nodes(name=pa)[0].get_descendants():
                # print "yay"
                for j in range(dist - 1):
                    # print M[g], S.search_nodes(name=x)[0].name, dist
                    # print dist
                    T[(g_child, S.search_nodes(name=x)[0].up.name)] = ['L', (g_child, x), (None, None)]
                    x = S.search_nodes(name=x)[0].up.name
                g_value.append((g_child, x))
            else:
                # print dist
                for j in range(dist - 1):
                    # print pa, S.search_nodes(name=x)[0].name, dist
                    # print dist
                    next_node = [node.name for node in S.search_nodes(name=x)[0].children \
                        if S.search_nodes(name=pa)[0] in S.search_nodes(name=node.name)[0].get_descendants() \
                        or pa == node.name]
                        
                    # print [node.name for node in S.search_nodes(name=x)[0].children]
                    # print next_node
                    T[(g_child, next_node[0])] = ['L', (g_child, x), (None, None)]
                    x = next_node[0]
                g_value.append((g_child, x))
                
        if len(g_value) == 1:
            g_value += [(None, None), (None, None)]
        T[(g, M[g])] = g_value
        
    return T
                

def pull_up_gene_node(G, S, alpha, g):
    L, M, sigma, delta, theta, xi, tau = alpha
    
    if g in sigma:
        sigma.remove(g)
        delta.append(g)
        # print 'sigma'
    elif g in delta:
        M[g] = S.search_nodes(name=M[g])[0].up.name
        # print 'delta'
    elif g in theta:
        # print 'theta'
        M[g] = S.search_nodes(name=M[g])[0].up.name
        if S.get_common_ancestor(M[g], tau[g]).name == M[g] and M[g] != tau[g]:
            # print "theta to delta"
            theta.remove(g)
            delta.append(g)
            g1, g2 = G.search_nodes(name=g)[0].children
            if (g, g1.name) in xi:
                xi.remove((g, g1.name))
            else:
                xi.remove((g, g2.name))
            del tau[g]
            
    return alpha
        

def explore(F, visited, g, source):
    visited[g] = True
    for node in F[g]:
        if node is not None:  # note that a leaf in F has children list [None], not []
            if node == source:
                return True
            if node not in visited:
                if explore(F, visited, node, source):
                    return True
    return False
    
    
def is_cycle(F, g):
    return explore(F, {}, g, g)
    
    
def find_first_cycle(G, G_dict, S, S_dict, alpha, flag=0):
    L, M, sigma, delta, theta, xi, tau = alpha
    T = dtl_to_recon_tree(S, G, alpha)
    F = buildReconciliation(S_dict, G_dict, T)
    # F_reformatted = nx.from_dict_of_lists(F).to_directed()
    
    post_order_S = [node.name for node in S.traverse("postorder")]
    nodes = sorted([node.name for node in G.search_nodes() if node not in G.search_nodes(children=[])], \
        key=lambda g: post_order_S.index(M[g]))
    # nodes = sorted([node.name for node in G.search_nodes() if node not in G.search_nodes(children=[])], \
    #     key=lambda g: 0 if S.get_tree_root().name == M[g] else S.get_distance(S.get_tree_root().name, M[g]), reverse=True)
    # cycles = nx.simple_cycles(F_reformatted)
    if flag == 1:
        print 'F', F
        print
        print post_order_S
        print
        print [M[g] for g in nodes]
    # for g in nodes:
    #     if g in [node for cycle in cycles for node in cycle]:
    #         return g
    for g in nodes:
        if is_cycle(F, g):
            return g
            
    return None
    
    
def temporal_consistency_fixer(G, G_dict, S, S_dict, alpha):
    g = find_first_cycle(G, G_dict, S, S_dict, alpha)
    i = 0
    while g is not None:
        i += 1
        # print i, g, alpha[1][g]
        alpha = pull_up_gene_node(G, S, alpha, g)
        # print 'new', g, alpha[1][g]
        g = find_first_cycle(G, G_dict, S, S_dict, alpha)
    return alpha, i
        

def out(S, G, alpha, outFile):
    # print alpha
    # print
    T = dtl_to_recon_tree(S, G, alpha)
    # print T
    # print
    d, s, t, l = 0, 0, 0, 0
    for key in T.keys():
        if T[key][0] == 'D':
            d += 1
        elif T[key][0] == 'S':
            s += 1
        elif T[key][0] == 'T':
            t += 1
        elif T[key][0] == 'L':
            l += 1
    print "D:", d, "S:", s, "T:", t, "L:", l, "total:", d * 2 + t * 3 + l
    outFile.write("D: " + str(d) +  " S: " + str(s) + " T: " + str(t) + " L: " + str(l) + \
        " total: " + str(d * 2 + t * 3 + l) + "\n")

    return d * 2 + t * 3 + l


def preprocess(st):

    while st.find(":") != -1 :
        i = st.find(":")
        j = i + 1
        while st[j] not in [";", ",", ")"] :
            j += 1
        st = st[: i] + st[j: ]
    return st


def eteTreeReader(fileName):
    """ Takes a fileName as input and returns the hostTree and parasiteTree in ETE Tree format"""
    
    fileHandle = open(fileName, 'r')
    contents = fileHandle.read()
    fileHandle.close()

    hostString, parasiteString, phiString = contents.split(";")
    hostString = hostString.strip()
    parasiteString = parasiteString.strip()
    hostString += ";"
    parasiteString += ";"
    # print hostString
    # print
    # print parasiteString
    hostString = preprocess(hostString)
    parasiteString = preprocess(parasiteString)

    # print hostString
    # print
    # print parasiteString

    hostTree = Tree(hostString, format=8)
    parasiteTree = Tree(parasiteString, format=8)

    return hostTree, parasiteTree


def main() :
    dVal = 2
    tVal = 3
    lVal = 1

    for i in xrange(100):

        index = str(i + 1)
        for j in xrange(4 - len(str(i + 1))):
            index = "0" + index
    
        fileName = "real-100taxa/COG" + index + ".newick"
        if not os.path.isfile(fileName):
            continue

        outFile = open("outputs/COG" + index + ".txt", 'w')

        print fileName[13:]
        outFile.write(fileName[13:] + "\n")

        S_dict, G_dict, _ = newickFormatReader.getInput(fileName)
        S, G = eteTreeReader(fileName)
        recs, allRecs = MasterReconciliation.Reconcile(["", fileName, str(dVal), str(tVal), str(lVal), "unit", "0", "1", "0", "1"])

        totRecs = len(allRecs)
        # print S
        # print G
        # print 
        # print S_dict
        # print
        # print G_dict
        # print

        print "# of Infeasible Reconciliations:", len(recs)
        outFile.write("# of Reconciliations: " + str(totRecs) + "\n")
        outFile.write("# of Infeasible Reconciliations: " + str(len(recs)) + "\n")
        
        min_cost = None

        for T in recs:
            # print T
            alpha = recon_tree_to_dtl(T)
            out(S, G, alpha, outFile)
            alpha, pull_up = temporal_consistency_fixer(G, G_dict, S, S_dict, alpha)
            cost = out(S, G, alpha, outFile)
            if min_cost is None or cost < min_cost:
                min_cost = cost
            print "number of operations:", pull_up
            outFile.write("number of operations: " + str(pull_up) + "\n")

        print "min total:", min_cost
        outFile.write("min total: " + str(min_cost) + "\n")

        outFile.close()


if __name__ == "__main__" :
   main()
