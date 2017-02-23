(ns distributions.core
  (:import (org.apache.commons.math3.distribution AbstractIntegerDistribution AbstractRealDistribution IntegerDistribution RealDistribution NormalDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution GumbelDistribution LaplaceDistribution LevyDistribution LogisticDistribution LogNormalDistribution NakagamiDistribution NormalDistribution ParetoDistribution TDistribution TriangularDistribution UniformRealDistribution WeibullDistribution BinomialDistribution EnumeratedIntegerDistribution GeometricDistribution HypergeometricDistribution PascalDistribution PoissonDistribution UniformIntegerDistribution ZipfDistribution MultivariateNormalDistribution)))

(defn beta
  [alpha beta]
  (new BetaDistribution alpha beta))
(defn binomial
  [n p]
  (new BinomialDistribution n p))
(defn cauchy
  [location scale]
  (new CauchyDistribution location scale))
(defn chi-squared
  [df]
  (new ChiSquaredDistribution df))
(defn degenerate
  [x]
  (new ConstantRealDistribution x))
(defn discrete-real
  [x p]
  (let [x-double (double-array x)
        p-double (double-array p)]
    (new EnumeratedRealDistribution x-double p-double)))
(defn discrete-integer
  [x p]
  (let [x-int (int-array x)
        p-double (double-array p)]
    (new EnumeratedIntegerDistribution x-int p-double)))
(defn exponential
  [mean]
  (new ExponentialDistribution mean))
(defn f-distribution
  [df1 df2]
  (new FDistribution df1 df2))
(defn gamma
  [shape scale]
  (new GammaDistribution shape scale))
(defn geometric
  [p]
  (new GeometricDistribution p))
(defn gumbel
  [location scale]
  (new GumbelDistribution location scale))
(defn hypergeometric
  [population-size number-successes sample-size]
  (new hypergeometric population-size number-successes sample-size))
(defn laplace
  [location scale]
  (new LaplaceDistribution location scale))
(defn double-exponential
  [location scale]
  (laplace location scale))
(defn levy
  [location scale]
  (new LevyDistribution location scale))
(defn logistic
  [location scale]
  (new LogisticDistribution location scale))
(defn log-normal
  [scale shape]
  (new LogNormalDistribution scale shape))
                                        ;(defn mvnormal)
(defn nakagami
  [shape spread]
  (new NakagamiDistribution shape spread))
(defn negative-binomial
  [r p]
  (new PascalDistribution r p))
(defn normal
  [location scale]
  (new NormalDistribution location scale))
(defn pareto
  [scale shape]
  (new ParetoDistribution scale shape))
(defn pascal
  [r p]
  (negative-binomial r p))
(defn poisson
  [mean]
  (new PoissonDistribution mean))
(defn t-distribution
  [df]
  (new TDistribution df))
(defn triangular
  [a c b]
  (new TriangularDistribution a c b))
(defn uniform-integer
  [lower upper]
  (new UniformIntegerDistribution lower upper))
(defn uniform-real
  [lower upper]
  (new UniformRealDistribution lower upper))
(defn weibull
  [shape scale]
  (new WeibullDistribution shape scale))
(defn zipf
  [number-of-elements exponent]
  (new ZipfDistribution number-of-elements exponent))


(defprotocol univariate-real
  (cdf [d] [d x])
  (pdf [d] [d x])
  (log-pdf [d] [d x])
  (mean [d])
  (variance [d])
  (support-lower-bound [d])
  (support-upper-bound [d])
  (icdf [d] [d x])
  (support-connected [d])
  (sample [d] [d n]))

(extend-type RealDistribution
  univariate-real
  (cdf
    ([d x] (.cumulativeProbability d x))
    ([d] (fn [x] (cdf d x))))
  (pdf
    ([d x] (.density d x))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (.logDensity d x))
    ([d] (fn [x] (log-pdf d x))))
  (mean [d] (.getNumericalMean d))
  (variance [d] (.getNumericalVariance d))
  (support-lower-bound [d] (.getSupportLowerBound d))
  (support-upper-bound [d] (.getSupportUpperBound d))
  (icdf
    ([d x] (.inverseCumulativeProbability d x))
    ([d] (fn [x] (icdf d x))))
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  )

(extend-type AbstractRealDistribution
  univariate-real
  (cdf
    ([d x] (.cumulativeProbability d x))
    ([d] (fn [x] (cdf d x))))
  (pdf
    ([d x] (.density d x))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (.logDensity d x))
    ([d] (fn [x] (log-pdf d x))))
  (mean [d] (.getNumericalMean d))
  (variance [d] (.getNumericalVariance d))
  (support-lower-bound [d] (.getSupportLowerBound d))
  (support-upper-bound [d] (.getSupportUpperBound d))
  (icdf
    ([d x] (.inverseCumulativeProbability d x))
    ([d] (fn [x] (icdf d x))))
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  )

(defprotocol univariate-integer
  (cdf [d] [d x])
  (pmf [d] [d x])
  (log-pmf [d] [d x])
  (mean [d])
  (variance [d])
  (support-lower-bound [d])
  (support-upper-bound [d])
  (icdf [d] [d x])
  (support-connected [d])
  (sample [d] [d n]))

(extend-type IntegerDistribution
  univariate-integer
  (cdf
    ([d x] (.cumulativeProbability d x))
    ([d] (fn [x] (cdf d x))))
  (pmf
    ([d x] (if (integer? x)
             (.probability d x)
             0))
    ([d] (fn [x] (pdf d x))))
  (log-pmf
    ([d x] (if (integer? x)
             (.logProbability d x)
             Double/NEGATIVE_INFINITY))
    ([d] (fn [x] (log-pdf d x))))
  (mean [d] (.getNumericalMean d))
  (variance [d] (.getNumericalVariance d))
  (support-lower-bound [d] (.getSupportLowerBound d))
  (support-upper-bound [d] (.getSupportUpperBound d))
  (icdf
    ([d x] (.inverseCumulativeProbability d x))
    ([d] (fn [x] (icdf d x))))
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  )

(extend-type AbstractIntegerDistribution
  univariate-integer
  (cdf
    ([d x] (.cumulativeProbability d x))
    ([d] (fn [x] (cdf d x))))
  (pmf
    ([d x] (if (integer? x)
             (.probability d x)
             0))
    ([d] (fn [x] (pdf d x))))
  (log-pmf
    ([d x] (if (integer? x)
             (.logProbability d x)
             Double/NEGATIVE_INFINITY))
    ([d] (fn [x] (log-pdf d x))))
  (mean [d] (.getNumericalMean d))
  (variance [d] (.getNumericalVariance d))
  (support-lower-bound [d] (.getSupportLowerBound d))
  (support-upper-bound [d] (.getSupportUpperBound d))
  (icdf
    ([d x] (.inverseCumulativeProbability d x))
    ([d] (fn [x] (icdf d x))))
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  )

(defprotocol exponential-family
  (natural-parameters [dist] )
  (play [dist] (.getMean dist)))
