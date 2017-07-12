For univariate distributions which extend protocols implementing cdf and icdf, generating independent draws can be achieved using the inverse-cdf method. For those extending only the protocol implementing cdf, icdf can fall back onto its default root finding algorithm, and the inverse-cdf method can still be used.

For other more exotic distributions of interest other sampling algorithms may be needed.

To illustrate these algorithms lets assume a trivial distribution that has a density propotional to x^2 on the unit interval

```
λ> (defn target-density [x]
     (cond
       (> x 1) 0
       (< x 0) 0
       :else (/ (* x x) 3)))
```





