# Posterior Predictive
The posterior predictive is implemented as the composition of posterior and marginalization operations. 

Here are some examples...

## Normal Normal
```
λ> (def data (sample (normal 3 1) 100))
#'λ/data

λ> (posterior-predictive data (normal :mu 1) (normal 0 5))
#distributions.core.Normal{:mean 2.833833314407506, :variance 506/501}
```

## Poisson Gamma
```
λ> (def data (sample (poisson 10) 100))
#'λ/data

λ> (posterior-predictive data (poisson :rate) (gamma 2 0.1))
#distributions.core.NegativeBinomial{:failures 962, :probability 0.990108803165183}
```

## Dirichlet Process
```
λ> (posterior-predictive [1.1 2.2 2.2] :G (dirichlet-process 10 (normal 0 1)))
#distributions.core.Mixture{:components [#distributions.core.DiscreteReal{:locations (1.1 2.2), :probabilities (1/3 2/3)} #distributions.core.Normal{:mean 0, :variance 1}], :probabilities (1/11 10/11)}
```
