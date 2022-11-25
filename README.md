---
title: "NMF R package -- Algorithms and Framework for Nonnegative Matrix Factorization (NMF)"
---

<!-- badges: start -->
  [![R-CMD-check](https://github.com/renozao/NMF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/renozao/NMF/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Background

Nonnegative Matrix Factorization (NMF) is an unsupervised learning technique that has been applied successfully in several fields, including signal processing, face recognition and text mining.
Recent applications of NMF in bioinformatics have demonstrated its ability to extract meaningful information from high-dimensional data such as gene expression microarrays. Developments in NMF theory and applications have resulted in a variety of algorithms and methods.
However, most NMF implementations have been on commercial platforms, while those that are freely available typically require programming skills.
This limits their use by the wider research community.

## Results
Our objective is to provide the bioinformatics community with an open-source, easy-to-use and unified interface to standard NMF algorithms, as well as with a simple framework to help implement and test new NMF methods.
For that purpose, we have developed a package for the R/BioConductor platform. The package ports public code to R, and is structured to enable users to easily modify and/or add algorithms.
It includes a number of published NMF algorithms and initialization methods and facilitates the combination of these to produce new NMF strategies.
Commonly used benchmark data and visualization methods are provided to help in the comparison and interpretation of the results.

## Conclusions
The NMF package helps realize the potential of Nonnegative Matrix Factorization, especially in bioinformatics, providing easy access to methods that have already yielded new insights in many applications.

## Availability

Documentation, source code and sample data are available in different places:

* Latest stable release from CRAN: https://cran.r-project.org/package=NMF

* Development versions on Github: 
	* Stable: https://github.com/renozao/NMF@master
	* Bleed-edge: https://github.com/renozao/NMF@devel

-----
