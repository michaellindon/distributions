# Integration and Expectation

## Expectation
`(expectation d)`

`(expectation d & options)`

Options are `:n`, a positive integer (default 10000) and `:method`, either "monte-carlo" or "quantile-integration" (default "monte-carlo").

Taking expectations of functions w.r.t. probability measure is decoupled into two parts for performance reasons.

```
λ> (expectation (normal 3 2))
#function[distributions.core/expectation-mc/fn--35747]

λ> (def E_N32 (expectation (normal 3 2)))
#'λ/E_N32

λ> (E_N32 identity)
2.9840497981514345

λ> (E_N23 (fn [x] (* (- x 3) (- x 3))))
1.9547371857236546
```
First note that `(expectation (normal 3 2))` returns a function. This can then be later applied to functions of which we wish to take the expectation w.r.t. a `(normal 3 2)` distribution.

The `expectation` function accepts two optional keyword arguments, namely, `n` and `method` which default to 10000 and "monte-carlo". An alternative to monte carlo is "quantile-integration" which is more accurate than monte carlo, but limited to 1 dimensional distributions. The trade-off between speed and accuracy can be controlled by the keyword argument `n`. Internally when `expectation` is called, a random sample ("monte-carlo") or grid ("quantile-integration") of size `n` is generated and stored so that it need not be recomputed when calling the resulting function multiple times.

```
λ> (def E_N32 (expectation (normal 3 2) :method "quantile-integration"))
#'λ/E_N32

λ> (E_N32 identity)
3.0010404833273636

λ> (E_N32 (fn [x] (* (- x 3) (- x 3))))
2.0073645213141385
```

```
λ> (def E_N32 (expectation (normal 3 2) :method "quantile-integration" :n 100000))
#'λ/E_N32

λ> (E_N32 identity)
3.0000981920167717

λ> (E_N32 (fn [x] (* (- x 3) (- x 3))))
2.0005300751367643
```

## Kullback-Leibler Divergence
Computing the KL-divergence between two probability distributions is implemented similarly

```
λ> (kullback-leibler (normal 3 2))
#function[distributions.core/kullback-leibler/fn--35763]

λ> ((kullback-leibler (normal 3 2)) (normal 4 2))
0.24797805571777043

λ> ((kullback-leibler (normal 3 2) :n 100 :method "quantile-integration") (normal 4 2))
0.2499999999999989
```
