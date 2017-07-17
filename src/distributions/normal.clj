(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.NormalDistribution)

(defrecord Normal [mean variance]
  random
  (sample [d] (.sample (new NormalDistribution mean (sqrt variance))))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new NormalDistribution mean (sqrt variance)) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (.cumulativeProbability (new NormalDistribution mean (sqrt variance)) x))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (.inverseCumulativeProbability (new NormalDistribution mean (sqrt variance)) x))
  (icdf [d] (fn [x] (icdf d x)))
  density-function
  (pdf [d x] (.density (new NormalDistribution mean (sqrt variance)) x))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d x] (.logDensity (new NormalDistribution mean (sqrt variance)) x))
  (log-pdf [d] (fn [x] (log-pdf d x)))
  support
  (support-lower [d] (.getSupportLowerBound (new NormalDistribution mean (sqrt variance))))
  (support-upper [d] (.getSupportUpperBound (new NormalDistribution mean (sqrt variance))))
  proximal
  (prox [d] (fn [h x] (* x (/ (inv h) (+ (inv h) (inv (variance d)))))))
  (prox [d h x] ((prox d) h x))
  first-moment
  (mean [d] (.getNumericalMean (new NormalDistribution mean (sqrt variance))))
  second-central-moment
  (variance [d] (.getNumericalVariance (new NormalDistribution mean (sqrt variance)))))

(defn normal [mean variance] (new Normal mean variance))

(defmethod posterior [distributions.core.Normal distributions.core.Normal]
  [data likelihood prior]
  (let [prior-var (:variance prior)
        prior-mean (:mean prior)
        obs-var (:variance likelihood)
        n (count data)
        post-var (inv (+ (/ n obs-var) (inv prior-var)))
        post-mean (* post-var (+
                               (/ prior-mean prior-var)
                               (* (reduce + data) (/ 1 obs-var))))]
    (normal post-mean post-var)))

(defmethod posterior [distributions.core.Normal distributions.core.InverseGamma]
  [data likelihood prior]
  (let [ss (reduce + 0 (map (fn [x] (* x x)) (map - data (repeat (:mean likelihood)))))
        n (count data)]
    (inverse-gamma (+ (:shape prior) (/ n 2)) (+ (:scale prior) (/ ss 2)))))

(defmethod marginal [distributions.core.Normal distributions.core.Normal]
  [likelihood prior]
  (normal (:mean prior) (+ (:variance likelihood) (:variance prior))))
