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
  (new HypergeometricDistribution population-size number-successes sample-size))
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

(defprotocol probability-function
  (P [d] [d x]))
(extend-protocol probability-function
  IntegerDistribution
  (P
    ([d] (fn [x] (.probability d x)))
    ([d x] (.probability d x)))
  RealDistribution
  (P
    ([d] (fn [x] (.probability d x)))
    ([d x] (.probability d x)))
  )

(defprotocol density-function
  (pdf [d] [d x])
  (log-pdf [d] [d x]))
(extend-protocol density-function
  RealDistribution
  (pdf
    ([d x] (.density d x))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (.logDensity d x))
    ([d] (fn [x] (log-pdf d x)))))

(defprotocol mass-function
  (pmf [d] [d x])
  (log-pmf [d] [d x]))
(extend-protocol mass-function
  IntegerDistribution
  (pmf
    ([d x] (if (integer? x)
             (.probability d x)
             0))
    ([d] (fn [x] (pdf d x))))
  (log-pmf
    ([d x] (if (integer? x)
             (.logProbability d x)
             Double/NEGATIVE_INFINITY))
    ([d] (fn [x] (log-pmf d x))))
  )

(defprotocol distribution-function
  (cdf [d] [d x]))
(extend-protocol distribution-function
  IntegerDistribution
  (cdf
    ([d] (fn [x] (.cumulativeProbability d x)))
    ([d x] (.cumulativeProbability d x)))
  RealDistribution
  (cdf
    ([d] (fn [x] (.cumulativeProbability d x)))
    ([d x] (.cumulativeProbability d x)))
  )

(defprotocol inverse-distribution-function
  (icdf [d] [d x]))
(extend-protocol inverse-distribution-function
  IntegerDistribution
  (icdf
    ([d] (fn [x] (.inverseCumulativeProbability d x)))
    ([d x] (.inverseCumulativeProbability d x)))
  RealDistribution
  (icdf
    ([d] (fn [x] (.inverseCumulativeProbability d x)))
    ([d x] (.inverseCumulativeProbability d x)))
  )

(defprotocol support
  (support [d])
  (support-lower [d])
  (support-upper [d])
  )
(extend-protocol support
  IntegerDistribution
  (support [d] [(support-lower d) (support-upper d)])
  (support-lower [d] (.getSupportLowerBound d))
  (support-upper [d] (.getSupportUpperBound d))
  RealDistribution
  (support [d] [(support-lower d) (support-upper d)])
  (support-lower [d] (.getSupportLowerBound d))
  (support-upper [d] (.getSupportUpperBound d))
  )

(defprotocol first-moment
  (mean [d]))
(extend-protocol first-moment
  IntegerDistribution
  (mean [d] (.getNumericalMean d))
  RealDistribution
  (mean [d] (.getNumericalMean d))
  )

(defprotocol second-central-moment
  (variance [d]))
(extend-protocol second-central-moment
  IntegerDistribution
  (variance [d] (.getNumericalVariance d))
  RealDistribution
  (variance [d] (.getNumericalVariance d))
  )

(defprotocol random
  (sample [d] [d n]))
(extend-protocol random
  IntegerDistribution
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  RealDistribution
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  )



(defprotocol exponential-family
  (natural-parameters [dist] )
  (play [dist] (.getMean dist)))

(defn arithmetic-mean [coll]
  (if (empty? coll)
    nil
    (let [sum (reduce + 0 coll)
          n (count coll)]
      (/ sum n))))


(Math/exp 3)
(defn quantile-integrate
  ([f d n]
   (let [delta (/ 1.0 n)
         grid (map (icdf d) (range delta 1 delta))
         f_i (map f grid)]
     (/ (reduce + 0 f_i) n))))

(sample (normal 0 1))

(def N (normal 0 1))
(def sample-size 10000)
(def draws (sort (sample N sample-size)))
(def cdfs (map (cdf N) draws))
(def zips (map vector (drop 1 cdfs) cdfs))
(def deltas (map (fn [x] (- (first x) (second x))) zips))

(sample (normal 0 1))
