# Posterior Distributions
`(posterior data likelihood prior)`

The `posterior` function can be used to get the posterior distribution in cases where the prior is conjugate to the likelihood. A keyword is used to express upon which parameter the prior is placed.

Here are some examples...

## Normal Normal
```
λ> (def data (sample (normal 3 1) 100))
#'λ/data

λ> (posterior data (normal :mu 1) (normal 0 5))
#distributions.core.Normal{:mean 3.120829760129205, :variance 5/501}
```
## Poisson Gamma
```
λ> (def data (sample (poisson 10) 100))
#'λ/data

λ> (posterior data (poisson :rate) (gamma 2 0.1))
#distributions.core.Gamma{:shape 992, :rate 100.1}
```

## Dirichlet Process
```
λ> (posterior [1.1 2.2 2.2] :G (dirichlet-process 10 (normal 0 1)))
#distributions.core.DirichletProcess{:concentration 13, :base-measure #distributions.core.Mixture{:components [#distributions.core.DiscreteReal{:locations (1.1 2.2), :probabilities (1/3 2/3)} #distributions.core.Normal{:mean 0, :variance 1}], :probabilities (1/11 10/11)}}
```

## TODO
The `posterior` function is implemented as a Clojure multimethod which dispatches on the type of the likelihood, the type of the prior, and in which likelihood parameter(s) keyword(s) appear.
The dispatching function is currently not sophisticated enough to dispatch on mathematical expressions such as `(normal (+ a (* c :mu)) 3)`, yet this should be possible using [expresso](https://github.com/clojure-numerics/expresso) and is on the TODO list.
