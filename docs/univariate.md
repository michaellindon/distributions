# Univariate Distributions
The following should implement `sample`, `pdf`, `log-pdf`, `support-lower`, `support-upper`, `cdf`, `icdf`, `mean`, `variance` and `pmf` and `log-pmf` for discrete distributions.

## Dirichlet
`(dirichlet concentration)`
```
λ> (dirichlet [1 10 20])
#distributions.core.DirichletDistribution{:alpha [1 10 20]}
```
## Discrete (Integer)
Expect better performance than using `discrete-real` with integers.
`(discrete-integer integers probabilities)`
```

λ> (discrete-integer [3 4] [0.5 0.5])
#distributions.core.DiscreteInteger{:integers [3 4], :probabilities (0.5 0.5)}
```

## Discrete (Real)
`(discrete-real locations probabilities)`
```
λ> (discrete-real [0.2 0.3] [0.5 0.5])
#distributions.core.DiscreteReal{:locations [0.2 0.3], :probabilities (0.5 0.5)}
λ> 
```

## Enumerated
Generalization of both Discrete-Real and Discrete-Integer
`(enumerated items probabilities)`
```
λ> (enumerated [:a "foo" 5 [1 2] {:mykey "mymap"}] [0.2 0.2 0.2 0.2 0.2])
#distributions.core.Enumerated{:items [:a "foo" 5 [1 2] {:mykey "mymap"}], :probabilities (0.2 0.2 0.2 0.2 0.2), :discrete #distributions.core.DiscreteInteger{:integers (0 1 2 3 4), :probabilities (0.2 0.2 0.2 0.2 0.2)}, :imap {:a 0, "foo" 1, 5 2, [1 2] 3, {:mykey "mymap"} 4}}
λ> (sample (enumerated [:a "foo" 5 [1 2] {:mykey "mymap"}] [0.2 0.2 0.2 0.2 0.2]))
{:mykey "mymap"}
```
## Gamma
`(gamma shape rate)`
```
λ> (gamma 3 2)
#distributions.core.Gamma{:shape 3, :rate 2}
```
## Generalized Double Pareto
`(gdp scale shape)`
```
λ> (gdp 4 3)
#distributions.core.GeneralizedDoublePareto{:scale 4, :shape 3}
```
## Inverse Gaussian
`(inverse-gaussian mean shape)`
```
λ> (inverse-gaussian 3 4)
#distributions.core.InverseGaussianDistribution{:mean 3, :shape 4}
```
## Negative Binomial
`(negative-binomial failues probabilities)`
```
λ> (negative-binomial 10 0.4)
#distributions.core.NegativeBinomial{:failures 10, :probability 0.4}
```
## Normal-Laplace
The density of this distribution is proportional to the product of a normal density with mean mu and variance gamma*s2 and a laplace density with rate tau. It is implemented as a mixture of truncated distributions. Useful in Bayesian computations involving normal observations and Laplace priors.

`(normal-laplace mu s2 gamma tau)`
```
λ> (normal-laplace 0 1 1 1)
#distributions.core.Mixture{:components [#distributions.core.TruncatedDistribution{:distribution #distributions.core.Normal{:mean 1.0, :variance 1.0}, :lower -Infinity, :upper 0, :F-lower 0.0, :F-upper 0.15865525393145702} #distributions.core.TruncatedDistribution{:distribution #distributions.core.Normal{:mean -1.0, :variance 1.0}, :lower 0, :upper Infinity, :F-lower 0.841344746068543, :F-upper 1.0}], :probabilities [0.5 0.5]}
```

## Normal/Gaussian
`(normal mean variance)`
```
λ> (normal 3 2)
#distributions.core.Normal{:mean 3, :variance 2}
```
## Poisson
`(poisson rate)`
```
λ> (poisson 1)
#distributions.core.Poisson{:rate 1}
```
## Polya-Gamma
`(polya-gamma z)`
```
λ> (polya-gamma 3)
#distributions.core.PolyaGammaDistribution{:z 3}
```
## Students t 
Standard:
`(t-distribution df)`

Location-Scale:
`(t-distribution df location scale)`
```
λ> (t-distribution 1 2 3)
#distributions.core.LocationScaleDistribution{:distribution #object[org.apache.commons.math3.distribution.TDistribution 0x38b1a404 "org.apache.commons.math3.distribution.TDistribution@38b1a404"], :location 2, :scale 3}
```
