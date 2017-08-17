(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.GammaDistribution)

(defrecord InverseGamma [shape scale]
  random
  (sample [d] (/ 1 (.sample (new GammaDistribution shape (/ 1 scale)))))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new GammaDistribution shape (/ 1 scale)) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (if (<= x 0) 0.0  (- 1 (.cumulativeProbability (new GammaDistribution shape (/ 1 scale)) (/ 1 x)))))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (if (= x 0) 0.0 (/ 1 (.inverseCumulativeProbability (new GammaDistribution shape (/ 1 scale)) (- 1 x)))))
  (icdf [d] (fn [x] (icdf d x)))
  density-function
  (pdf [d x] (if (<= x 0) 0.0 (* (/ 1 (* x x)) (.density (new GammaDistribution shape (/ 1 scale)) (/ 1 x)))))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d x] (if (<= x 0) Double/NEGATIVE_INFINITY (+ (* -2 (Math/log x)) (.logDensity (new GammaDistribution shape (/ 1 scale)) (/ 1 x)))))
  (log-pdf [d] (fn [x] (log-pdf d x)))
  support
  (support-lower [d] 0.0)
  (support-upper [d] Double/POSITIVE_INFINITY)
  first-moment
  (mean [d] (/ scale (dec shape)))
  second-central-moment
  (variance [d] (/ (* scale scale) (* (dec shape) (dec shape) (- shape 2)))))

(defn inverse-gamma [shape scale] (new InverseGamma shape scale))
