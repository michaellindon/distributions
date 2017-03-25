(ns distributions.core
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.linear :as la])
  (:import (org.apache.commons.math3.distribution AbstractIntegerDistribution AbstractRealDistribution IntegerDistribution RealDistribution NormalDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution GumbelDistribution LaplaceDistribution LevyDistribution LogisticDistribution LogNormalDistribution NakagamiDistribution NormalDistribution ParetoDistribution TDistribution TriangularDistribution UniformRealDistribution WeibullDistribution BinomialDistribution EnumeratedIntegerDistribution GeometricDistribution HypergeometricDistribution PascalDistribution PoissonDistribution UniformIntegerDistribution ZipfDistribution)))

(set-current-implementation :vectorz)

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
  [shape rate]
  (new GammaDistribution shape (/ 1 rate)))
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
(defn logistic-distribution
  [location scale]
  (new LogisticDistribution location scale))
(defn log-normal
  [scale shape]
  (new LogNormalDistribution scale shape))
(defn nakagami
  [shape spread]
  (new NakagamiDistribution shape spread))
(defn negative-binomial
  [r p]
  (new PascalDistribution r p))
(defrecord MultivariateNormalDistribution [mean variance])
(defn mvnormal
  ([mean variance]
   (MultivariateNormalDistribution. mean variance))
  ([variance]
   (let [n (column-count variance)
         mean (zero-vector n)]
     (mvnormal mean variance))))
(defn normal
  [mean variance]
  (new NormalDistribution mean (Math/sqrt variance)))
(defn pareto
  [scale shape]
  (new ParetoDistribution scale shape))
(defn pascal
  [r p]
  (negative-binomial r p))
(defn poisson
  [mean]
  (new PoissonDistribution mean))
(defrecord LSTDistribution [location scale df])
(defn t-distribution
  ([df]
   (new TDistribution df))
  ([location scale df]
   (LSTDistribution. location scale df)))
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
    ([d] (fn [x] (log-pdf d x))))
  LSTDistribution
  (pdf
    ([d x] (/ (pdf (new TDistribution (:df d)) (/ (- x (:location d)) (:scale d))) (:scale d)))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (- (log-pdf (new TDistribution (:df d)) (/ (- x (:location d)) (:scale d))) (Math/log (:scale d)) ))
    ([d] (fn [x] (log-pdf d x))))
  MultivariateNormalDistribution
  (log-pdf
    ([d x]
     (let [mean-vector (:mean d)
           cov-matrix (:variance d)
           n (ecount x)
           residual (sub x mean-vector)
           exponent (* 0.5 (dot residual (mmul (inverse cov-matrix) residual)))
           normalizer (+ (* 0.5 (log (det cov-matrix)))
                         (* n 0.5 (log (* 2 Math/PI))))]
       (negate (+ normalizer exponent))
       ))
    ([d] (fn [x] (log-pdf d x))))
  (pdf
    ([d x] (Math/exp (log-pdf d x)))
    ([d] (fn [x] (pdf d x))))
  )

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
  LSTDistribution
  (icdf
    ([d x] (+ (:location d) (* (:scale d) (.inverseCumulativeProbability (new TDistribution (:df d)) x))))
    ([d] (fn [x] (icdf d x)))
    )
  )

(defprotocol support
  (support-interval [d])
  (support-lower [d])
  (support-upper [d])
  )
(extend-protocol support
  IntegerDistribution
  (support-interval [d] [(support-lower d) (support-upper d)])
  (support-lower [d] (.getSupportLowerBound d))
  (support-upper [d] (.getSupportUpperBound d))
  RealDistribution
  (support-interval [d] [(support-lower d) (support-upper d)])
  (support-lower [d] (.getSupportLowerBound d))
  (support-upper [d] (.getSupportUpperBound d))
  )

(defprotocol first-moment
  (mean [d]))
(extend-protocol first-moment
  clojure.lang.Sequential
  (mean [coll] (if (empty? coll)
              nil
              (let [sum (reduce + coll)
                    n (count coll)]
                (/ sum n))))
  IntegerDistribution
  (mean [d] (.getNumericalMean d))
  RealDistribution
  (mean [d] (.getNumericalMean d))
  MultivariateNormalDistribution
  (mean [d] (:mean d))
  )

(defn square [x] (* x x))
(defn negate [x] (* -1 x))
(def probit (cdf (normal 0 1)))

(defprotocol second-central-moment
  (variance [d]))
(extend-protocol second-central-moment
  clojure.lang.Sequential
  (variance [coll] (if (empty? coll)
                 nil
                 (let [mu (mean coll)]
                   (mean (map (fn [x] (square (- x mu))) coll)))))
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
  LSTDistribution
  (sample
    ([d] (+ (:location d) (* (:scale d) (sample (new TDistribution (:df d))))))
    ([d n] (take n (repeatedly #(sample d)))))
  MultivariateNormalDistribution
  (sample
    ([d]
     (let [cov-matrix (:variance d)
           mean-vector (:mean d)
           ncols (column-count cov-matrix)
           z (sample (normal 0 1) ncols)
           d (la/svd cov-matrix)
           {U :U S :S V* :V*} (la/svd cov-matrix)
           D (diagonal-matrix (map (fn [x] (sqrt (max x 0))) S))]
       (add mean-vector (mmul U D z))))
    ([d n] (take n (repeatedly #(sample d)))))
  )

(defn expectation
  ([h d n]
   (mean (map h (sample d n))))
  ([h d]
   (expectation h d 10000)))

(defn kullback-leibler
  ([f g n]
   (let [h (fn [x] (- (log-pdf f x) (log-pdf g x)))]
     (expectation h f n)))
  ([f g]
   (kullback-leibler f g 10000)))

(defn expectation-qi
  ([h d n]
   (let [delta (/ 1.0 n)
         grid (map (icdf d) (range delta 1 delta))]
     (mean (map h grid))))
  ([h d]
   (expectation-qi h d 10000)))

(defn kullback-leibler-qi
  ([f g n]
   (let [h (fn [x] (- (log-pdf f x) (log-pdf g x)))]
     (expectation-qi h f n)))
  ([f g]
   (kullback-leibler-qi 10000)))

(defn quantile-integrate
  ([f d n]
   (let [delta (/ 1.0 n)
         grid (map (icdf d) (range delta 1 delta))
         f_i (map f grid)]
     (/ (reduce + 0 f_i) n))))
