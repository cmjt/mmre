## Overiew

`R` package to fit two-stare continuous-time discrete-space Markov models with individual level random effects. 
Methodology developed with [Enrico Pirotta](mailto:pirotta.enrico@gmail.com) to assess the effects of exposure 
to Navy sonar on marine mammal movement patterns.

## Installation

To install the `mmre` package run

```r
devtools::install_github("cmjt/mmre")
```

The model likelihoods are coded using `TMB`, in oder to compile all `TMB` templates after installation run

```r
mmre::compile.mmre()
```

and to load the templates run 

```r
mmre::dll.mmre()
```

once in each workspace.