#!/usr/bin/env python
import sys
import random
from itertools import combinations
from Bio import SeqIO
import kmerizer
import numpy as np
import pydot



def fastaLoad(fileName):
    file = SeqIO.parse(open(fileName, 'r'), "fasta")
    scaffolds = [scaffold for scaffold in file]
    file.close()
    return scaffolds


def intersectionStats(nodePair, kmerSets, intersections):
    node1, node2 = nodePair
    k1 = kmerSets[node1]
    k2 = kmerSets[node2]
    pair = (node1, node2)
    if pair in intersections:
        stats = intersections[pair]
    else:
        #kint = np.intersect1d(k1, k2, assume_unique = True)
        common = kmerizer.count_common(k1, k2)
        stats = common
        intersections[pair] = stats
        #print pair, stats
        #sys.stdout.flush()
    return stats


def buildTree(tree, kmerSets, intersections):
    # calculate the intersection stats of kmers for all tree node pairs
    intStats = [(pair, intersectionStats(pair, kmerSets, intersections)) for pair in  combinations(tree, 2)]
    intStats.sort(lambda a, b: cmp(b[1], a[1]))
    closePair, statsx = intStats[0]
    n1, n2 = closePair
    k1 = kmerSets[n1]
    k2 = kmerSets[n2]
    kint = np.intersect1d(k1, k2, assume_unique = True)
    k1diff = np.setdiff1d(k1, k2, assume_unique = True)
    k2diff = np.setdiff1d(k2, k1, assume_unique = True)
    # build new tree: start by removing the nodes to be joined
    newTree = [x for x in tree if x not in (n1, n2)]
    newNode = (n1, n2)
    newTree.append(newNode)
    kmerSets[newNode] = kint
    kmerSets[n1] = k1diff
    kmerSets[n2] = k2diff
    #print n1, k1diff.size
    #print n2, k2diff.size
    print newNode, kint.size, k1diff.size, k2diff.size
    print 'tree:', newTree
    sys.stdout.flush()
    return newTree


def graphTree(graph, tree, kmerSets, parentNode):
    leaf = isinstance(tree, str)
    nodeLabel = str(kmerSets[tree].size)
    if leaf:
        nodeLabel += "\n" + tree
    node = pydot.Node(nodeLabel)
    graph.add_node(node)
    if parentNode:
        graph.add_edge(pydot.Edge(parentNode, node))
    if not leaf:
        graphTree(graph, tree[0], kmerSets, node)
        graphTree(graph, tree[1], kmerSets, node)

K = int(sys.argv[1])
hashSize = int(sys.argv[2])
suffix = '.genome.fa'
hashBits = random.getrandbits(K * 2)

kmerSets = {}
for arg in sys.argv[3:]:
    name = arg
    if name.endswith(suffix):
        name = name[:-len(suffix)]
    print 'Kmerizing', name
    sys.stdout.flush()
    scaffolds = fastaLoad(arg)
    kmerArray = np.unique(np.concatenate([kmerizer.kmerize(K, str(s.seq)) for s in scaffolds]))
    kmerArray ^= hashBits
    kmerArray.sort()
    # kmerSets[name] = kmerArray
    kmerSets[name] = kmerArray[:hashSize]

tree = kmerSets.keys()
nleaves = len(tree)
print 'Kmerized', nleaves, 'genomes'
sys.stdout.flush()

intersections = {}
while len(tree) > 1:
    tree = buildTree(tree, kmerSets, intersections)
    print len(tree)
    sys.stdout.flush()

graph = pydot.Dot(graph_type='graph', rankdir='LR')
graphTree(graph, tree[0], kmerSets, None)
graph.write_png("kmertree.png")
