(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.PoissonDistribution)

(defrecord Poisson [rate]
  random
  (sample [d] (.sample (new PoissonDistribution rate)))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new PoissonDistribution rate) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (.cumulativeProbability (new PoissonDistribution rate) x))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (.inverseCumulativeProbability (new PoissonDistribution rate) x))
  (icdf [d] (fn [x] (icdf d x)))
  mass-function
  (pmf [d] (fn [x] (pdf d x)))
  (pmf [d x] (if (integer? x) (.probability (new PoissonDistribution rate) x) 0))
  (log-pmf [d x] (if (integer? x) (.logProbability (new PoissonDistribution rate) x) Double/NEGATIVE_INFINITY))
  (log-pmf [d] (fn [x] (log-pmf d x)))
  density-function
  (pdf [d] (fn [x] (pdf d x)))
  (pdf [d x] (if (integer? x) (.probability (new PoissonDistribution rate) x) 0))
  (log-pdf [d x] (if (integer? x) (.logProbability (new PoissonDistribution rate) x) Double/NEGATIVE_INFINITY))
  (log-pdf [d] (fn [x] (log-pmf d x)))
  support
  (support-lower [d] (.getSupportLowerBound (new PoissonDistribution rate)))
  (support-upper [d] (.getSupportUpperBound (new PoissonDistribution rate)))
  first-moment
  (mean [d] (.getNumericalMean (new PoissonDistribution rate)))
  second-central-moment
  (variance [d] (.getNumericalVariance (new PoissonDistribution rate))))

(defn poisson [rate] (new Poisson rate))

(defmethod posterior [distributions.core.Poisson distributions.core.Gamma]
  [data likelihood prior]
  (gamma (+ (reduce + data) (:shape prior)) (+ (:rate prior) (count data))))

(defmethod marginal [distributions.core.Poisson distributions.core.Gamma]
  [likelihood {b1 :shape b2 :rate }]
  (negative-binomial b1 (/ b2 (inc b2))))
