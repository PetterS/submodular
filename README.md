[![Build Status](https://travis-ci.org/PetterS/submodular.png)](https://travis-ci.org/PetterS/submodular)

Please refer to Readme.pdf for documentation. Automatically generated documentation is available in doc/html.

This library implements the optimization method decribed in (1). The code may be freely used as long as this publication is cited.

* (1) Fredrik Kahl and Petter Strandmark, ["Generalized Roof Duality for Pesudo-Boolean Optimization"](http://www2.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=482), International Conference on Computer Vision, 2011.
* (2) Fredrik Kahl and Petter Strandmark, ["Generalized Roof Duality"](http://www2.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=520), Discrete Applied Mathematics, 2012.
* (3) Petter Strandmark, ["Discrete Optimization in Early Vision"](http://www2.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=537), Phd thesis.

It also contains an implementation of 

* (4) Alexander Fix, Aritanan Gruber, Endre Boros, Ramin Zabih, "A Graph Cut Algorithm for Higher-order Markov Random Fields", International Conference on Computer Vision, 2011.

and it is able to use

* (5) Hiroshi Ishikawa, "Transformation of General Binary MRF Minimization to the First Order Case", PAMI 2011.
* (6) C. Rother, V. Kolmogorov, V. Lempitsky, and M. Szummer, "Optimizing binary MRFs via extended roof duality", CVPR 2007.

if the software is downloaded. If you use the methods/software from (4–6), you should cite them.

Building
--------
The library is built using CMake. Users of Ubuntu can look at the .travis.yml file, which contains all commands neccessary for downloading all requirements and building the library.

On Windows, the library builds with later versions of Visual Studio.

The current build status on Ubuntu is:
[![Build Status](https://travis-ci.org/PetterS/submodular.png)](https://travis-ci.org/PetterS/submodular)

Reproducing figures 1 and 2 in the paper
----------------------------------------
 * Compile the demo program
 * Run:
   * batchrun3
   * batchrun4
 * From MATLAB, run:
   * plot_batchrun('run_3_1000_1000.data')
   * plot_batchrun('run_4_1000_300.data')

This will reproduce the figures from (2) (similar to (1)).

 * (optional) Use Python to run:
   * python parse_batchrun.py

 This will generate statistics.


Running the examples from the papers
------------------------------------
 The file "examples from paper.cmd" will run the two examples
 which appear in the papers.
 
![Screenshot](http://www.strandmark.net/submodular.png)
