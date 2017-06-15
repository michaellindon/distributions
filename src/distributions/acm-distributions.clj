(in-ns 'distributions.core)

(import '(org.apache.commons.math3.distribution AbstractIntegerDistribution AbstractRealDistribution IntegerDistribution RealDistribution NormalDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution GumbelDistribution LaplaceDistribution LevyDistribution LogisticDistribution LogNormalDistribution NakagamiDistribution NormalDistribution ParetoDistribution TDistribution TriangularDistribution UniformRealDistribution WeibullDistribution BinomialDistribution EnumeratedIntegerDistribution GeometricDistribution HypergeometricDistribution PascalDistribution PoissonDistribution UniformIntegerDistribution ZipfDistribution))

(defn beta
  [alpha beta]
  (new BetaDistribution alpha beta))
(defn binomial
  [n p]
  (new BinomialDistribution n p))
(defn bernoulli
  [p]
  (binomial 1 p))
(defn cauchy
  [location scale]
  (new CauchyDistribution location scale))
(defn chi-squared
  [df]
  (new ChiSquaredDistribution df))
(defn degenerate
  [x]
  (new ConstantRealDistribution x))
(defn discrete-integer
  [x p]
  (let [x-int (int-array x)
        p-double (double-array p)]
    (new EnumeratedIntegerDistribution x-int p-double)))
(defn exponential
  [rate]
  (new ExponentialDistribution (/ 1 rate)))
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
  [location rate]
  (new LaplaceDistribution location (/ 1 rate)))
(defn double-exponential
  [location rate]
  (laplace location rate))
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
(defn normal
  [mean variance]
  (new NormalDistribution mean (sqrt variance)))
(defn pareto
  [scale shape]
  (new ParetoDistribution scale shape))
(defn pascal
  [r p]
  (negative-binomial r p))
(defn poisson
  [mean]
  (new PoissonDistribution mean))
(defn triangular
  [a c b]
  (new TriangularDistribution a c b))
(defn uniform-integer
  [lower upper]
  (new UniformIntegerDistribution lower upper))
(defn uniform
  [lower upper]
  (new UniformRealDistribution lower upper))
(defn weibull
  [shape scale]
  (new WeibullDistribution shape scale))
(defn zipf
  [number-of-elements exponent]
  (new ZipfDistribution number-of-elements exponent))

(defn get-rate [d] (/ 1 (.getScale d)))
(defn get-scale [d] (.getScale d))

(extend-protocol probability-function
  IntegerDistribution
  (P
    ([d] (fn [x] (.probability d x)))
    ([d x] (.probability d x)))
  RealDistribution
  (P
    ([d] (fn [x] (.probability d x)))
    ([d x] (.probability d x))))

(extend-protocol density-function
  RealDistribution
  (pdf
    ([d x] (.density d x))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (.logDensity d x))
    ([d] (fn [x] (log-pdf d x)))))

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
    ([d] (fn [x] (log-pmf d x)))))

(extend-protocol distribution-function
  IntegerDistribution
  (cdf
    ([d] (fn [x] (.cumulativeProbability d x)))
    ([d x] (.cumulativeProbability d x)))
  RealDistribution
  (cdf
    ([d] (fn [x] (.cumulativeProbability d x)))
    ([d x] (.cumulativeProbability d x))))

(extend-protocol support
  IntegerDistribution
  (support-lower [d] (.getSupportLowerBound d))
  (support-upper [d] (.getSupportUpperBound d))
  RealDistribution
  (support-lower [d] (.getSupportLowerBound d))
  (support-upper [d] (.getSupportUpperBound d)))

(extend-protocol inverse-distribution-function
  IntegerDistribution
  (icdf
    ([d] (fn [x] (.inverseCumulativeProbability d x)))
    ([d x] (.inverseCumulativeProbability d x)))
  RealDistribution
  (icdf
    ([d] (fn [x] (.inverseCumulativeProbability d x)))
    ([d x] (.inverseCumulativeProbability d x))))

(extend-protocol first-moment
  IntegerDistribution
  (mean [d] (.getNumericalMean d))
  RealDistribution
  (mean [d] (.getNumericalMean d)))

(extend-protocol second-central-moment
  IntegerDistribution
  (variance [d] (.getNumericalVariance d))
  RealDistribution
  (variance [d] (.getNumericalVariance d)))


(extend-protocol proximal
  NormalDistribution
  (prox
    ([d] (fn [h x] (* x (/ (inv h) (+ (inv h) (inv (variance d)))))))
    ([d h x] ((prox d) h x)))
  LaplaceDistribution
  (prox
    ([d] (fn [h x] (* (signum x) (max 0.0 (- (abs x) (* h (get-rate d)))))))
    ([d h x] ((prox d) h x))))

(extend-protocol random
  IntegerDistribution
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d)))))
  RealDistribution
  (sample
    ([d] (.sample d))
    ([d n] (take n (repeatedly #(sample d))))))

