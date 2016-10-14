#!/usr/bin/env python
import sys
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
        kint = np.intersect1d(k1, k2, assume_unique = True)
        stats = (kint.size, float(kint.size) / float(min(k1.size, k2.size)))
        intersections[pair] = stats
    return stats


def buildTree(tree, kmerSets, intersections):
    # calculate the intersection stats of kmers for all tree node pairs
    intStats = [(pair, intersectionStats(pair, kmerSets, intersections)) for pair in  combinations(tree, 2)]
    intStats.sort(lambda a, b: cmp(b[1][0], a[1][0]))
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
    print n1, k1diff.size
    print n2, k2diff.size
    print newNode, kint.size
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
suffix = '.genome.fa'

kmerSets = {}
for arg in sys.argv[2:]:
    name = arg
    if name.endswith(suffix):
        name = name[:-len(suffix)]
    print 'Kmerizing', name
    sys.stdout.flush()
    scaffolds = fastaLoad(arg)
    kmerArray = np.unique(np.concatenate([kmerizer.kmerize(K, str(s.seq)) for s in scaffolds]))
    kmerSets[name] = kmerArray

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
