
<!-- README.md is generated from README.Rmd. Please edit that file -->

# {Rprobit}: Estimation of Probit Models for discrete choice data (ordered and unordered choices) including the panel case.

<!-- badges: start -->
<!-- badges: end -->

## Disclaimer 

The package is in early stages of development. It may contain bugs and 
broken functionality. The software comes with no warranties whatsoever. 

## Introduction 
{Rprobit} implements the specification and estimation of multinomial
probit models that are used in the context of discrete choice data sets.
Thereby individuals (named deciders) decide between a finite number of
distinct alternatives. Both the alternatives as well as the deciders may
be characterized using some regressor variables that influence the
choices.

The package includes code for the cross-sectional case where each
decider takes only one choice as well as the panel case wherein deciders
face a number of choice occasions. The number of choice occasions may
vary among deciders as can the number of alternatives presented to them.

Multinomial probit models (MNP) represent one version of the random
utility type models that assume that deciders take their decisions based
on an underlying utility assigned to each alternative. The alternative
with the highest utility is then chosen. The utilities are modelled as
linear functions of the regressor variables. Individual heterogeneities
in preferences can be modelled using the concept if mixing of
parameters. {Rprobit} implements Gaussian mixing of the parameters.

Beside the classical MNP models for unordered choice data wherein each
alternative is assigned a separate utility function {Rprobt} can also
handle ordered choice situations wherein the alternatives are ordered
such as the answers to survey questions or the amount of pieces bought
from a product. In this case only one utility function per choice is
modelled.

## Installation

You can install the released version of {RprobitB} from the ZIP file.

## Documentation

The package is documented in several vignettes.

## Example

We analyze a data set of 2929 stated choices by 235 Dutch individuals
deciding between two virtual train trip options based on the price, the
travel time, the level of comfort, and the number of changes. The data
is saved in the {mlogit} package. We transform the travel time from
minutes to hours and the travel price from guilders to euros:

``` r
data("Train", package = "mlogit")
Train$price_A <- Train$price_A / 100 * 2.20371
Train$price_B <- Train$price_B / 100 * 2.20371
Train$time_A <- Train$time_A / 60
Train$time_B <- Train$time_B / 60
str(Train)
#> 'data.frame':    2929 obs. of  11 variables:
#>  $ id       : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ choiceid : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ choice   : Factor w/ 2 levels "A","B": 1 1 1 2 2 2 2 2 1 1 ...
#>  $ price_A  : num  52.9 52.9 52.9 88.1 52.9 ...
#>  $ time_A   : num  2.5 2.5 1.92 2.17 2.5 ...
#>  $ change_A : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ comfort_A: num  1 1 1 1 1 0 1 1 0 1 ...
#>  $ price_B  : num  88.1 70.5 88.1 70.5 70.5 ...
#>  $ time_B   : num  2.5 2.17 1.92 2.5 2.5 ...
#>  $ change_B : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ comfort_B: num  1 1 0 0 0 0 1 0 1 0 ...
```

The following lines fit a probit model that explains the chosen trip
alternatives (`choice`) by their `price`, `time`, number of `change`s,
and level of `comfort` (the lower this value the higher the comfort).
For normalization, the first linear coefficient, the `price`, is fixed
to `-1`, which allows to interpret the other coefficients as monetary
values:

``` r
# estimate a MNL model using mlogit
Train.ml <- mlogit:::mlogit.data(Train, choice = "choice", shape = "wide",
                     varying = 4:11, alt.levels = c("A", "B"), sep = "_",
                     opposite = c("price", "time", "change", "comfort"))
#> Warning in dfidx::dfidx(data = data, dfa$idx, drop.index = dfa$drop.index, : the
#> levels shouldn't be provided with a data set in wide format

Train.mod <- mlogit:::mlogit(choice ~ price + time + change + comfort | 0, Train.ml, reflevel = "A", R=500)

summary(Train.mod)
#> 
#> Call:
#> mlogit:::mlogit(formula = choice ~ price + time + change + comfort | 
#>     0, data = Train.ml, reflevel = "A", R = 500, method = "nr")
#> 
#> Frequencies of alternatives:choice
#>       A       B 
#> 0.50324 0.49676 
#> 
#> nr method
#> 5 iterations, 0h:0m:0s 
#> g'(-H)^-1g = 0.00014 
#> successive function values within tolerance limits 
#> 
#> Coefficients :
#>          Estimate Std. Error z-value  Pr(>|z|)    
#> price   0.0673580  0.0033933 19.8506 < 2.2e-16 ***
#> time    1.7205514  0.1603517 10.7299 < 2.2e-16 ***
#> change  0.3263409  0.0594892  5.4857 4.118e-08 ***
#> comfort 0.9457256  0.0649455 14.5618 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-Likelihood: -1724.2
# set up the probit model
train_mod <- setup_Rprobit(form = choice ~ price + time + change + comfort | 0 ,  data_raw = Train,  ids = c("id"))

# change scale fixation to price parameter. 
#train_mod$mod$lthL = 1
train_mod$mod$HL = matrix(0,3,1)
train_mod$mod$HL[3,1] = 1
train_mod$mod$fL = matrix(0,3,1)

#train_mod$mod$lthb = 3
train_mod$mod$Hb = diag(4)[,-1]
train_mod$mod$fb = as.matrix(c(-0.06,0,0,0),ncol=1)

# initialise with MNL results
train_mod$theta_0 <- -c(as.numeric(Train.mod$coefficients[2:4]),1)

train_fit <- fit_Rprobit(Rprobit_obj = train_mod,init_method = "theta", 
                         control_nlm = list(approx_method = "SJ", probit = FALSE, hess = 0))
```

The estimated effects can be presented via:

``` r
summary(train_fit)
#> Summary of probit model:
#> 
#> N = 235 decision makers
#> alt = 2 alternatives
#> Tp = 5 to 19 choice occasions
#> 
#> U = X * beta + eps
#> beta = b
#> eps ~ MVN(0,Sigma)
#> 
#> b =
#> price    -0.060000 (0.000000) 
#> time     -1.550696 (0.126086) 
#> change   -0.295149 (0.051991) 
#> comfort  -0.866770 (0.052271) 
#> 
#> Sigma =
#>    A                    B                    
#> A  0.000000 (0.000000)  0.000000 (0.000000) 
#> B  0.000000 (0.000000)  2.332485 (0.194281) 
#> 
#>                                 A         B
#> true choice percentages 0.5032434 0.4967566
#> averages of predictions 0.4966261 0.5033739
#> 
#> Log-likelihood: -1727.694945
```
