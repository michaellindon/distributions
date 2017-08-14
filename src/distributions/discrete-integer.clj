(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution)

(defrecord DiscreteInteger [integers probabilities]
  random
  (sample [d] (.sample (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities))))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (.cumulativeProbability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (.inverseCumulativeProbability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x))
  (icdf [d] (fn [x] (icdf d x)))
  mass-function
  (pmf [d] (fn [x] (pdf d x)))
  (pmf [d x] (if (integer? x) (.probability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x) 0))
  (log-pmf [d x] (if (integer? x) (.logProbability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x) Double/NEGATIVE_INFINITY))
  (log-pmf [d] (fn [x] (log-pmf d x)))
  density-function
  (pdf [d] (fn [x] (pdf d x)))
  (pdf [d x] (if (integer? x) (.probability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x) 0))
  (log-pdf [d x] (if (integer? x) (.logProbability (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)) x) Double/NEGATIVE_INFINITY))
  (log-pdf [d] (fn [x] (log-pmf d x)))
  support
  (support-lower [d] (.getSupportLowerBound (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities))))
  (support-upper [d] (.getSupportUpperBound (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities))))
  first-moment
  (mean [d] (.getNumericalMean (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities))))
  second-central-moment
  (variance [d] (.getNumericalVariance (new EnumeratedIntegerDistribution (int-array integers) (double-array probabilities)))))

(defn discrete-integer
  [integers probabilities & {:keys [log?]
                             :or {log? false}}]
  (new DiscreteInteger integers (if log? (normalize-log probabilities) (normalize probabilities))))
