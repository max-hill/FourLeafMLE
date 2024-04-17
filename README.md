---
abstract: |
  We solve a global optimization problem to perform statistical
  inference using DNA sequence alignment data for quartets under the
  Cavender-Farris-Neyman model from phylogenetics.
author:
- |
  Max Hill\
  Department of Mathematics\
  University of California, Riverside\
  Riverside, USA, 92521\
  [max.hill1@ucr.edu](max.hill1@ucr.edu){.uri}
- |
  Jose Israel Rodriguez\
  Department of Mathematics\
  Department of Electrical and Computer Engineering\
  University of Wisconsin --- Madison\
  Madison, USA, 53703\
  [jose@math.wisc.edu](jose@math.wisc.edu){.uri}
bibliography:
- refs.bib
title: A maximum likelihood estimator for quartets under the
  Cavender-Farris-Neyman model
titlehead: A maximum likelihood estimator for quartets under the CFN
  model
---

# Introduction: tree topologies and the tree of life

One of the core motivations in evolutionary biology is to reconstruct
the *tree of life* that describes how species have evolved over time. A
common approach is to use DNA sequence alignment data from present-day
taxa to infer a *phylogenetic tree*, understood as a leaf-labelled tree
with edges representing ancestral populations, vertices representing
common ancestors, and leaves representing extant taxa. The edge weights,
or *branch lengths* represent some measure of evolutionary time (usually
expected number of mutations per site), and the tree *topology*
describes clade structure of the taxa. For instance,
Figure [1](#fig:long-branch-attraction){reference-type="ref"
reference="fig:long-branch-attraction"} has four leaves that are labeled
in such a way to indicate species $1$ and $2$ have a more recent common
ancestor with one another than with species $3$ or $4$.

A phylogenetic tree thus represents a hypothesis about the evolutionary
history of a set of taxa; the statistical inference problem is to infer
the tree topology and branch lengths from a DNA sequence of fixed
length. Our implementation solves polynomial systems to perform
this inference.

# Maximum likelihood estimation and quartet methods

A standard approach to inferring a phylogenetic tree is to use maximum
likelihood estimation. By summarizing data as a vector $(u_1,\dots,u_n)$
of counts, the goal is to maximum the likelihood of the data. In other
words, we seek the value of $\theta$ that maximizes
$p_1^{u_1}\cdots p_n^{u_n}$ where $p_i$ is the probability of observing
event $i$ and each $p_i$ depends on the parameters $\theta$ of the
model. The maximum likelihood problem addressed in `fourLeafMLE`
involves estimating the branch lengths and unrooted tree topology of a
4-leaf tree from sequence data of a fixed length $k$ generated according
to the Cavendar-Farris-Neyman (CFN) model [@semple2003phylogenetics].

Considerable research has focused on understanding the properties of
maximum likelihood estimation. Even in the simplest cases of 3- and
4-leaf trees, the problem exhibits substantial complexity, with a
general solution known for 3-leaf trees, but not 4-leaf trees
[@allman2019maximum; @chor2007analytic; @CS2004; @GGS2022; @hill2024maximum; @hobolth2024maximum; @KK2019; @Y2000].
The $4$-leaf case considered here is of special interest first due to
the popularity of quartet-based inference methods for inferring both
phylogenetic trees and networks [@warnow2017book], and second because it
is the simplest case in which the phenomenon of *long-branch attraction*
can be observed, a form of estimation bias which is only partially
understood [@B2005; @SR2021].

Two challenges in optimizing the likelihood function are its non-convex
nature and the presence of numerous boundary cases during optimization.
Moreover, there can be multiple distinct maximizers of the likelihood
function, and the maximizers can have branch lengths which are
infinitely long or zero [@CHHP2000; @steel-1994]. This package resolves
these two obstacles for the CFN model by using computer algebra and the
theory of maximum likelihood (ML) degrees [@HKS2005]. One feature
distinguishing this software is that it implements a framework to
characterize case when the tree estimates have infinite or zero-length
branches. We believe this is useful in obtaining an understanding of
ways that maximum likelihood inference can fail, and as a first step
leads to Figure [3](#fig:the-figure){reference-type="ref"
reference="fig:the-figure"}.

# Running the maximum likelihood estimator

Our main contribution puts algebra into practice by implementing a
custom tailored solver for algebraic statisticians. Given data in the
form of a site frequency vector, the main functionality is to return a
list of global maximizers of the likelihood function. The output
includes information about each of the maximizers, for example the tree
configuration and optimal branch lengths.

::: leftbar
    SITE_PATTERN_DATA = [212, 107, 98, 115, 114, 89, 102, 163]
    fourLeafMLE(SITE_PATTERN_DATA) 
    ## Output
    1-element Vector{Vector{Any}}:
     [-2036.1212979788797, 
      "R1", 
      [1], 
      [0.1380874841, 0.46347429951, 0.5231552324, 0.3975875363, 0.6835395124], 
      "θ1, θ2, θ3, θ4, θ5", 
      "binary quartet with topology τ=1"]
:::

For complete details on the code, see the code documentation provided here, but now we
give a broad overview about our software. The CFN model is a
parameterized semialgebraic set in $\mathbb{R}^8$ with dimension $6$. We
optimize the likelihood function by partitioning the semialgebraic set
according to boundaries from the parameterization, similiar to the
approach in [@allman2019maximum] for a different model. There are $10$
different boundaries up to symmetry. All but one of them have ML degree
less than two. The last one has ML degree 92 while the main component of
the semialgebraic set has ML degree fourteen.[^1] We optimize the
likelihood function on each of the boundary cases by solving the
likelihood equations [@HKS2005] by using analytic expressions when the
ML degree is less than two and
`Homotopy Continuation.jl` [@HomotopyContinuation.jl] to compute the
critical points otherwise.

# Putting into practice: Visualizing the geometry of data

Our method is fast. On a desktop machine with 12 i5-10400 CPUs
(2.90GHz), the `fourLeafMLE` solves on average 50 global optimization
problems per second. This is sufficiently fast to create high-quality
visualizations of the geometry of the maximum likelihood problem. For
example, Figure
[2](#fig:topology-classification-hadamard-plot){reference-type="ref"
reference="fig:topology-classification-hadamard-plot"} presents a grid
of $90,000$ points, each representing a choice of Hadamard edge length
parameters $\theta_x,\theta_y$ (see [@semple2003phylogenetics]) for the
tree shown in Figure
[1](#fig:long-branch-attraction){reference-type="ref"
reference="fig:long-branch-attraction"}. For each point, maximum
likelihood estimation was performed using random data (1000b) drawn from
the model, and the point colored according to which binary quartet
topologies were found to be compatible with the maximum likelihood
estimates.

<figure id="fig:the-figure">
<figure id="fig:long-branch-attraction">

<figcaption>A quartet with topology <span
class="math inline">12|34</span>. When <span
class="math inline"><em>x</em></span> is very small and <span
class="math inline"><em>y</em></span> is large, the maximum likelihood
estimate tends to have topology <span class="math inline">23|14</span>,
a bias is known as long-branch attraction.</figcaption>
</figure>
<figure id="fig:topology-classification-hadamard-plot">
<img
src="images/classification-plot-hadamard-k1000-resolution300x300-marker0.57.png" />
<figcaption>MLE topologies for random data, represented by
color.</figcaption>
</figure>
<figcaption>A visualization of 90,000 inferred tree topologies (right)
for random data consisting of 1000bp drawn according to the CFN model on
trees of the form given on the left. The reddish bottom right corner of
(b) provides a new visualization of the long-branch attraction
phenomenon.</figcaption>
</figure>

[^1]: The numbers $14$ and $92$ were computed previously in [@GGS2022
    Table 2].
