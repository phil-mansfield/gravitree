# gravitree
A tree code for quickly computing gravitational forces for large collections of particles. Written in Go. Documentation can be found [here](https://pkg.go.dev/github.com/phil-mansfield/gravitree).

## Background

Computing gravitational forces typically requires O(N^2) calculations between every pair of particles. There are many approaches to getting these forces without doing this much work. One popular approach is to build a tree which breaks particles into a hierarchy groups based on position. Examples of these types of trees include [octrees](https://en.wikipedia.org/wiki/Octree), [k-d trees](https://en.wikipedia.org/wiki/K-d_tree), [R trees](https://en.wikipedia.org/wiki/R-tree), and uniformly spaced bins.

Once the tree is built, the user can compute the gravitational force/potential at a point by exploring the tree. Starting with the root node, the user checks if the particles inside that node are suffiicently far away that they can be represented by a single particle. (This is called a "monopole". You can also apply more accurate "multipole" approximations.) If so, the gravitational force is calculated from the center of mass of that cell. If not, the node is "opened" and the user plays the same game with the child nodes. The decision to open a node is made using an "opening criteria."

The original/classic version of this approach is the [Barnes-Hut](https://ui.adsabs.harvard.edu/abs/1986Natur.324..446B/abstract) algorithm. This algorithm uses an octree as the underlying data structure and uses an "opening angle" as its opening criteria (i.e. the ratio between node size and the distance to the node). Extremely advanced codes are based on this approach (The most notable of which is probably [Springel et al. 2001](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract)). There are many alternative approaches to contructing this tree that have been introduced over the years.

There there are currently no gravitational tree codes wirtten in Go that I am aware of.

## Technical Description

`gravitree` is designed to accurately compute gravitational forces in individual bound objects containing ~10^3 - 10^7 particles, such as dark matter haloes in simulations.

`gravitree` uses a k-d tree as the core data structure, similar to [PKDGRAV](https://arxiv.org/abs/1609.08621). k-d trees are simple to implement, fast to construct, and tend to have a slimmer memory footprint than octtrees due to the smaller number of empty/sparsely populated nodes. k-d trees can have large depths and some constructions can lead to oblong node shapes which have large errors in their monopole moments. These problems are mitigated by other design choices in `gravitree`. 

`gravitree` splits nodes halfway along the longest dimension until the leaf size of 16 points is reached. Other approaches (such as taking the median or center of mass position in the node) lead to more oblong node shapes and thus makes multipole approximations innaccurate at larger distances. Closed nodes are appoximated as monopoles and nodes are opened using the PKDGRAV3 opening criteria. To reduce tree traversal time, opening decisions are made collectively for each leaf node in the tree. 

All these choices (other than the underlying data structure) can be altered at runtime.

## Performance

`gravitree` is a fairly fast code and outperforms many similar C and Fortran codes, even when configured similarly. I've optimized `gravitree`'s hot loops a decent amount, and `gravitree` uses a more cache- and allocation-friendly memory format than many of its peers. For example, when configured identically to the `Rockstar` halo finder's `fast3tree`, `gravitree` produces potentials for Einasto point distributions in about 60% the time. Additionally, the default parameters have been extensively tested to minimize runtime while maintaining high force accuracy.

However, `gravitree` does not use GPU acceleration and does not use multi-node parallelism and thus cannot approach truly massive problems. It's best used for smaller problems, such as halo finding or small single-timestep simulations.

I haven't added thread paralellization inside the `gravitree` routines yet (although this would not be difficult to do), but given that this code is designed to work on ~< 10^7 particle objects, the user can just set up a worker queue and run each tree as a task in the queue. This approach would essentially lead to prefect parallelism for large object groups.
