# Geting Started

## Creating and Getting Parameters of Distributions
Lets start with an example. To create a Normal distribution with mean 1 and variance 2 we write
```
λ> (normal 1 2)
#distributions.core.Normal{:mean 1, :variance 2}
```
Distributions are implemented as records. As seen above, the `normal` function called with arguments 1 and 2 returns a Normal record with mean parameter 1 and variance parameter 2. Parameters of the distribution can be accessed by using the appropriate keywords.

```
λ> (def mynormal (normal 0 1))
#'λ/mynormal

λ> mynormal
#distributions.core.Normal{:mean 0, :variance 1}

λ> (:mean mynormal)
0

λ> (:variance mynormal)
1

```
!!! note "Normal Parameterization" 
    The design choice has been made to parameterize the normal distribution by mean and variance instead of mean and stan dard deviation.

## Sampling Distributions
Sampling from a distribution is achieved using the `sample` function
```
λ> (sample mynormal)
0.26276240636941356

λ> (sample mynormal 3)
(0.04801945336086195 0.8026957377916983 0.22748199289423185)
```
Notice that `sample` has multiple arities. One can get a single draw using the 1-arity version, or many samples using the 2-arity version.
```
λ> (sample (poisson 4))
2

λ> (sample (poisson 4) 3)
(4 1 6)
```
The `sample` function dispatches against the distribution type (clojure record). This kind of polymorphism is achieved using Clojure protocols. Both arities of the `sample` function are implemented in the `random` protocol. A list of protocols can be found [here](https://github.com/michaellindon/distributions/blob/master/src/distributions/protocols.clj). 

## Means and Variances
The Poisson is parameterized by its rate.
```
λ> (poisson 4)
#distributions.core.Poisson{:rate 4}

λ> (:rate (poisson 4))
4
```
Whilst the rate of a Poisson can be retrieved using the keyword `:rate`, other functions are required to get the mean and variance of this distribution. For this purpose there are the `mean` and `variance` functions.

```
λ> (mean (poisson 4))
4.0

λ> (variance (poisson 4))
4.0

λ> (standard-deviation (poisson 4))
2.0
```

These functions also work on `clojure.lang.Sequential` types.
```
λ> (mean '(1 5 3 2))
11/4

λ> (variance [9 3 1 5])
35/4

λ> (standard-deviation [4 3 1])
1.2472191289246473
```
## Density Evaluation
Probabilistic algorithms often require evaluations of a density or mass function. For this there are the functions `pdf`, `pmf`, `log-pdf` and `log-pmf`. The log-functions are not merely the pdf/pmf function composed with the logarithm function, but instead are computed on the log-scale. This is useful for avoiding underflow and overflow in computations. These functions have two arities. `(pdf d x)` evaluates the pdf of d at x, whereas `(pdf d)` returns a function, namely the pdf of d.
```
λ> (pdf (gamma 3 2) 1)
0.5413411329464509

λ> (log-pdf (gamma 3 2) 1)
-0.6137056388801094

λ> (map (pdf (gamma 3 2)) [1 2 3])
(0.5413411329464509 0.2930502222197469 0.08923507835998892)

λ> (map (log-pdf (gamma 3 2)) [1 2 3])
(-0.6137056388801094 -1.2274112777602189 -2.41648106154389)
```

!!! note "pdf or pmf":
    Whilst discrete distributions implement both pdf and pmf, continuous distributions only implement pdf. If in doubt just use `pdf` for all purposes. 
## Cumulative Distribution Function

The cumulative distribution function of a distribution d can be obtained using `(cdf d)`. Evaluating the cdf of a distribution d at x can be done directly with `(cdf d x)`.
```
λ> (cdf (exponential 1) 1)
0.6321205588285577

λ> (cdf (exponential 1))
#function[distributions.core/eval35150/fn--35151/fn--35152]

λ> ((cdf (exponential 1)) 1)
0.6321205588285577
```
The inverse or "quantile" function is used similarly.
```
λ> (icdf (exponential 1) 0.6321205588285577)
1.0

λ> (icdf (exponential 1))
#function[distributions.core/eval35174/fn--35175/fn--35176]
```

## Support
The support of the distribution can be obtained using `supprt-lower` and `support-upper`

```
λ> (support-lower (exponential 3))
0.0

λ> (support-upper (exponential 3))
Infinity
```

