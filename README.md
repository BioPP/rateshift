# Bio++ RateShift

A reimplementation of the RateShift test by Pupko and Galtier:
> A covarion-based method for detecting molecular adaptation: application to the evolution of primate mitochondrial genomes.
> Pupko T, Galtier N.
> Proc Biol Sci. 2002 Jul 7;269(1498):1313-6.

This implementation compares two sets of branches, arbitrarily refered to as 'foreground' and 'background'. The program compares for each site a model where the two sets evolve with a distinct evolutionary rate (estimated by the program) to a model where the two sets evolve under the same rate. The output file contains, for each site:
* the common rate
* the foreground and background rates
* the AIC of the one-rate and two-rate models
* the p-value of the likelihood ratio test between the two models

BppRateShift follows the BppO syntax, and all sequence and tree format supported by Bio++ (see the BppSuite manual for a complete reference). Input branches are specified as a list of node ids. Such nodes can be output by programms such as bppML, or visualized using the BppPhyView.

Installation
============

### Binaries

A binary static executable (linux 64bits) can be downladed [soon available].
    
### Compilation from sources

Bio++ Rateshift requires the latest development version of the Bio++ libraries (components core, seq and phyl). If they are not installed on your system, you can perform a local compilation:
```bash
git clone https://github.com/BioPP/bpp-core.git
cd bpp-core
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local
make install -j 2

git clone https://github.com/BioPP/bpp-seq.git
cd bpp-seq
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local
make install -j 2

git clone https://github.com/BioPP/bpp-phyl.git
cd bpp-phyl
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local
make install -j 2
```
If everything went well, you can compile Bio++ RateShift:
```bash
git clone https://github.com/BioPP/rateshift.git
cd rateshift
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local
make install
```

Example of usage
================

An example is provided in the `example` directory, using the *nd5* gene. An alignment with site selection (obtained with Gblocks) is in the `nd5.aln.mase` file. It can be visualized with the *seaview* software, or simply using a text editor. A phylogenetic tree generated with PhyML, using a LG08+Gamma(4) model is also provided and saved in the file `nd5-PhyML_tree.dnd`.

In order to select the foreground branches, we need to identify the ids of the corresponding branches. While bppRateShift only consider unrooted tree, a Newick representation is always rooted on a node, arbitrarily. Branches are then labelled internally by the id of the node down the branch:
```
    |-----0
----2
    |-----1
```
In this example, the two terminal branches have ids 0 and 1, and the branch leading to their common ancestor has id 2. It is possible to visualize the ids of each node by opening the tree in the Bio++ PhyView program, or to simply ask bppRateShift to output them in a standard Newick file:

```bash
bpprateshift param=rateshift.bpp output.tree_ids.file=nd5_ids.dnd
``` 
Here, the foreground branches are all the branches for the Primate clade, which correspond to nodes 13 to 23.

This information in hand, we can now run bppRateShift with the command:
```bash
bpprateshift param=rateshift.bpp foreground_branches=13-23
```

- Note 1: it is possible to specify a discontinuous range of branches, for instance `3,6,45,10` and to use a combined syntax like `3,5-8,34-41`.

- Note 2: bppRateShift can optimize branch lengths and model parameters prior to testing rate-shifts. This is recommended, as different software may have slightly different model parameterization. Optimizing model parameters will ensure that branch lengths are in the correct unit. This step may, however, be skipped by setting the `optimization=None` option. Optimization can be controled by the exact same parameters as the bppML program from the bppSuite software (but only homogeneous models are supported so far). Refer to the bppSuite manual for more details.

bppRateShift will test for all selected sites and output the results in tabular form, according to the file path specified in the output file.
