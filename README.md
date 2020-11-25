# gravitree
A tree code for quickly computing gravitational forces for large collections of particles. Written in Go.

Author: Phil Mansfield, Ph.D

## Background

Computing gravitational forces typically requires O(N^2) calculations between every pair of particles. There are many approaches to getting these forces without doing this much work with varying trade-offs. One popular approach is to build a tree which breaks particles into a hierarchy groups based on position. Examples of these types of trees include [octrees](https://en.wikipedia.org/wiki/Octree), [k-d trees](https://en.wikipedia.org/wiki/K-d_tree), [R trees](https://en.wikipedia.org/wiki/R-tree), and uniformly spaced bins.

Once the tree is built, the user can compute the gravitational force/potential at a point by exploring the tree. Starting with the root node, the user checks if the particles inside that node are suffiicently far away that they can be represented by a single particle. (This is called a "monopole". You can also apply more accurate "multipole" approximations.) If so, the gravitational force is calculated from the center of mass of that cell. If not, the node is "opened" and the user plays the same game with the child nodes. The decision to open a node is made using an "opening criteria."

The original/classic version of this approach is the [Barnes-Hut](https://ui.adsabs.harvard.edu/abs/1986Natur.324..446B/abstract) algorithm. This algorithm uses an octree as the underlying data structure and uses an "opening angle" as its opening criteria (i.e. the ratio between node size and the distance to the node). Extremely advanced codes are based on this approach (The most notable of which is probably [Springel et al. 2001](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract)). Later authors have pointed out that such criteria can lead to very large errors ([Salmon & Warren 1994](https://ui.adsabs.harvard.edu/abs/1994JCoPh.111..136S/abstract)).

Unfortunately, no such code exists in Go.

## Description

This code uses a k-d tree to partition space, monopole approximations in unopened nodes, and the maximum error relations in Salmon & Warren (1994) as opening criteria. At this level, the algorithm is identical to fast3tree from [Behroozi et al. (2013)](https://ui.adsabs.harvard.edu/abs/2013ApJ...762..109B/abstract). The implementation contains some of my own performance/memory optimizations as well as some optimizations from [Han et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..604H/abstract).
