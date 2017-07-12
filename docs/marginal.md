# Marginal Distributions

`(marginal likelihood prior)`

For compound distibutions, it is possible to compute the marginal distribution from the corresponding joint.
Suppose Y follows a (normal :mu 3) distribution and that mu follows a (normal 1 2) distribution.
A keyword to express which parameter is to be marginalized. 

Here are some examples...

## Normal Normal
```
λ> (marginal (normal :mu 3) (normal 1 2))
#distributions.core.Normal{:mean 1, :variance 5}
```
## Poisson Gamma
```
λ> (marginal (poisson :rate) (gamma 3 2))
#distributions.core.NegativeBinomial{:failures 3, :probability 2/3}
```

## Mixtures
This also works with mixture distributions
```
λ> (marginal (normal :mu 3) (mixture [(normal 1 2) (normal 3 4)] [0.1 0.9]))
#distributions.core.Mixture{:components (#distributions.core.Normal{:mean 1, :variance 5} #distributions.core.Normal{:mean 3, :variance 7}), :probabilities [0.1 0.9]}
```

## Dirichlet Process
A Dirichlet process is informally a distribution over distributions and so, adhering to the convention of expressing the object over which marginalization is to be carried out using a keyword, the marginal can be computed as

```
λ> (marginal :G (dirichlet-process 10 (normal 0 1)))
#distributions.core.Normal{:mean 0, :variance 1}
 
```

## TODO
The `marginal` function is implemented as a Clojure multimethod which dispatches on the type of the likelihood, the type of the prior, and in which likelihood parameter(s) keyword(s) appear.
The dispatching function is currently not sophisticated enough to dispatch on mathematical expressions such as `(normal (+ a (* c :mu)) 3)`, yet this should be possible using [expresso](https://github.com/clojure-numerics/expresso) and is on the TODO list.
