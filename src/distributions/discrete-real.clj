(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.EnumeratedRealDistribution)

(defrecord DiscreteReal [locations probabilities]
  random
  (sample [d] (.sample (new EnumeratedRealDistribution (double-array locations) (double-array probabilities))))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new EnumeratedRealDistribution (double-array locations) (double-array probabilities)) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (.cumulativeProbability (new EnumeratedRealDistribution (double-array locations) (double-array probabilities)) x))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (.inverseCumulativeProbability (new EnumeratedRealDistribution (double-array locations) (double-array probabilities)) x))
  (icdf [d] (fn [x] (icdf d x)))
  density-function
  (pdf [d x] (.density (new EnumeratedRealDistribution (double-array locations) (double-array probabilities)) x))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d x] (.logDensity (new EnumeratedRealDistribution (double-array locations) (double-array probabilities)) x))
  (log-pdf [d] (fn [x] (log-pdf d x)))
  support
  (support-lower [d] (.getSupportLowerBound (new EnumeratedRealDistribution (double-array locations) (double-array probabilities))))
  (support-upper [d] (.getSupportUpperBound (new EnumeratedRealDistribution (double-array locations) (double-array probabilities))))
  first-moment
  (mean [d] (.getNumericalMean (new EnumeratedRealDistribution (double-array locations) (double-array probabilities))))
  second-central-moment
  (variance [d] (.getNumericalVariance (new EnumeratedRealDistribution (double-array locations) (double-array probabilities)))))

(defn discrete-real [locations probabilities & {:keys [log?]
                                                :or {log? false}}] (new DiscreteReal
                                                                        locations
                                                                        (if log? (normalize-log probabilities) (normalize probabilities))))

(defmethod marginal [java.lang.Object distributions.core.DiscreteReal]
  [likelihood {locations :locations probabilities :probabilities}]
  (let [params (first (map first (filter (fn [[k v]] (keyword? v)) likelihood)))
        ]
    (mixture (map (fn [x] (assoc likelihood params x) ) locations) probabilities)))
