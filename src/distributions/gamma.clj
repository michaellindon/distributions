(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.GammaDistribution)

(defrecord Gamma [shape rate]
  random
  (sample [d] (.sample (new GammaDistribution shape (/ 1 rate))))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new GammaDistribution shape (/ 1 rate)) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (.cumulativeProbability (new GammaDistribution shape (/ 1 rate)) x))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (.inverseCumulativeProbability (new GammaDistribution shape (/ 1 rate)) x))
  (icdf [d] (fn [x] (icdf d x)))
  density-function
  (pdf [d x] (.density (new GammaDistribution shape (/ 1 rate)) x))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d x] (.logDensity (new GammaDistribution shape (/ 1 rate)) x))
  (log-pdf [d] (fn [x] (log-pdf d x)))
  support
  (support-lower [d] (.getSupportLowerBound (new GammaDistribution shape (/ 1 rate))))
  (support-upper [d] (.getSupportUpperBound (new GammaDistribution shape (/ 1 rate))))
  first-moment
  (mean [d] (.getNumericalMean (new GammaDistribution shape (/ 1 rate))))
  second-central-moment
  (variance [d] (.getNumericalVariance (new GammaDistribution shape (/ 1 rate)))))

(defn gamma [shape rate] (new Gamma shape rate))
