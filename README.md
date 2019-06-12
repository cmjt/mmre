
```
## Loading required package: TMB
```

```
## Loading required package: Rcpp
```

## Overiew

`R` package to fit two-stare continuous-time discrete-space Markov models with individual level random effects. 
Methodology was developed with [Enrico Pirotta](pirotta.enrico@gmail.com) to assess the effects of exposure 
to Navy sonar on marine mammal movement patterns.

## Installation

To install the `mmre` package run


```r
devtools::install_github("cmjt/mmre")
library(mmre)
```

All model likelihoods are coded using `TMB`; in oder to compile all `TMB` templates after installation run


```r
compile.mmre()
```

and to load the templates run 


```r
dll.mmre()
```

once in each workspace.

## Example

The `mmre` package contains an example dataset of estimated latitude and longitude locations 
of two individuals on the AUTEC Naval range. The states indicate if an individual was off 
(state = 1) or on (state = 2) range.

<img src="/figure/AUTEC.png" title="Estimated tracks on and around the AUTEC Naval range" alt="Estimated tracks on and around the AUTEC Naval range" width="30%" style="display: block; margin: auto;" />



```r
data(example)
data <- example
```
To fit a simple two state continuous-time Markov model run


```r
mod.basic <- fit.mmre(data = data,parameters = list(log_baseline = log(c(0.5,0.5))))
```

and to get the estimated transition probability matrix P(t = 1)


```r
get.probs(mod.basic,1)
>           State 1    State 2
> State 1 0.9326097 0.06739027
> State 2 0.4195685 0.58043148
```

To compare the results to the `msm` package run


```r
library(msm)
msm.fit <- msm(state ~ time, subject = ID, data = data, qmatrix = rbind(c(0, 0.5), c(0.5, 0)),  
    exacttimes = FALSE)
pmatrix.msm(msm.fit)
>           State 1    State 2
> State 1 0.9325607 0.06743932
> State 2 0.4194816 0.58051835
```

### Model with individual level random effects


```r
qs <- as.numeric(get.params(mod.basic,FALSE)[,1])
pars <- list(log_baseline = qs,u = matrix(numeric(2*length(table(data$ID))),
      ncol = 2),log_sigma = 0)
mod.basic.re <- fit.mmre(data = data,parameters = pars)
```


```r
get.probs(mod.basic.re,1)
>           State 1    State 2
> State 1 0.9333624 0.06663765
> State 2 0.5536407 0.44635932
```


```r
show.random(mod.basic.re)
```

<img src="figure/unnamed-chunk-12-1.png" title="Individual level random effects" alt="Individual level random effects" width="30%" />

