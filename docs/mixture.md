# Mixture Distributions
Mixture distributions can be created by specifying a vector of mixture components and a vector of probabilities

`(mixture components probabilities)`

```
λ> (mixture [(poisson 3) (gamma 2 5) (normal 0 4)] [0.1 0.3 0.6])
#distributions.core.Mixture{:components [#distributions.core.Poisson{:rate 3} #distributions.core.Gamma{:shape 2, :rate 5} #distributions.core.Normal{:mean 0, :variance 4}], :probabilities [0.1 0.3 0.6]}
```

The pdf at x is computed as a summation over products of component pdfs evaluated at x with component probabilities
```
(pdf (mixture [(poisson 3) (gamma 2 5) (normal 0 4)] [0.1 0.3 0.6]) 0.3)
0.6203866596063337
```

Because the pdf involves a sum over components, the log-pdf is computed using the log-sum-exp identity, so that the pdf at x for each component is still computed on a log scale.
```
λ> (log-pdf (mixture [(poisson 3) (gamma 2 5) (normal 0 4)] [0.1 0.3 0.6]) 0.3)
-0.47741235080208877
```

The cdf at x of a mixture is compted as a summation over products of component cdfs evaluated at x with component probabilities
```
λ> (cdf (mixture [(poisson 3) (gamma 2 5) (normal 0 4)] [0.1 0.3 0.6]) 0.3)
0.47340170214760935
```

There is no closed form expression for the inverse cdf of a mixture. Instead a root finding algorithm is used.

```
λ> (icdf (mixture [(poisson 3) (gamma 2 5) (normal 0 4)] [0.1 0.3 0.6]) 0.47340170214760935)
0.30000000000001137
```
