Tomé S. Silva 
====================
[ [Github](https://tomessilva.github.io/) ] [ [Gitlab](https://gitlab.com/t.s.silva/)  ] [ [ORCID](https://orcid.org/0000-0002-9434-8686) ] [ [Scholar](https://scholar.google.com/citations?hl=en&user=Z1dFRJEAAAAJ&view_op=list_works&sortby=pubdate) ]

* R packages

  * [fort](https://tomessilva.github.io/fort): Fast orthogonal random transforms in R. Provides convenient access to fast, structured, random linear transforms implemented in C++ that are (at least approximately) orthogonal or semi-orthogonal, and are often much faster than matrix multiplication. [Github repository](https://github.com/tomessilva/fort). [Package manual](https://tomessilva.github.io/manuals/fort_0.0.1.pdf).


* Numerical software

  * [femtogRad](https://gitlab.com/t.s.silva/femtograd): A minimal autograd library in R. Started as a port of Karpathy's micrograd Python library, with additional features being added on top. Supports training of simple neural networks for regression (least squares loss) and classification (max-margin loss) by backpropagation.
 
  * [GoDec](https://gitlab.com/t.s.silva/godec): R implementation of the GoDec algorithm, which decomposes a matrix (X) as a sum of a low-rank (L), a sparse (S) and a noise (N) matrix (X = L + S + N). Ported from the original Matlab code (GoDec.m, written by Tianyi Zhou) with small modifications.


* Web applications

  * [poweranalyser](https://sparos.shinyapps.io/poweranalyser_1/): A visual tool to perform power analyses for studies involving hypothesis testing and regression analysis.
  
  * [feedEst](https://webtools.sparos.pt/feedest/): A feeding rate prediction tool for gilthead seabream and European seabass, developed within the context of EU project [PerformFISH](https://web.archive.org/web/20230714003419/http://performfish.eu/).
  
  * [wastEst](https://webtools.sparos.pt/wastest/): wastEst is an application that provides reasonable estimates of waste outputs from fish rearing activities in tanks or cages, using minimal input from the user, to enable fair comparisons between scenarios in terms of their environmental impact. This tool was developed in the context of EU project [ARRAINA](http://web.archive.org/web/20210923211938/https://arraina.eu/).
  
  * [ficoEst](https://webtools.sparos.pt/ficoest/): A tool to estimate the body composition of farmed fish using incomplete measurements. Developed in the context of project [FICA](https://www.sparos.pt/projects/fica/). Supporting role (main developer: Filipe Soares).


* Proprietary software

  * [FEEDNETICS](https://www.sparos.pt/products/#feednetics): A commercial fish farm simulation tool based on a [state-of-the-art model of fish metabolism and growth](https://doi.org/10.3390/jmse11030472). Backend development.

  * fishR: An extensible framework for development, calibration and evaluation of ODE models of fish physiology and metabolism, implemented in R. Used internally at [SPAROS](https://www.sparos.pt/).
  
  * fishSim: A web application based on the fishR framework to enable calibration and simulation of fish growth models by end-users with no programming experience. Used internally at [SPAROS](https://www.sparos.pt/).
  
  * microfit: An R implementation of the [FiT algorithm](https://www.sparos.pt/products/#fit) for estimation of optimal feeding rates in fish farming, using the fishR framework. Used internally at [SPAROS](https://www.sparos.pt/).
  
  * bioact: Analysis framework to estimate *in silico* the bioactivity and sensorial properties of peptides present in feed ingredients. Used internally at [SPAROS](https://www.sparos.pt/).
  
  * sinkfeel: A web application to analyse the sinking profiles of extruded feed pellets. Used internally at [SPAROS](https://www.sparos.pt/).
  
  
* Other software and code

  * [2DExample](https://tomessilva.github.io/2DExample): Supplemental material for "[Visualization and differential analysis of protein expression data using R](http://web.archive.org/web/20220311094808/https://sapientia.ualg.pt/bitstream/10400.1/9779/1/manuscript_final.pdf)". Example proteomic dataset (from 2D gel electrophoresis) and R code.
