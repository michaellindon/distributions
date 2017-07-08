(in-ns distributions.core)

(defrecord drec [args])

(defn distribution [args] (drec. args))

(extend-protocol probability-function
  drec
  (P
    ([d] (fn [x] (.probability ACMDistribution x)))
    ([d x] (.probability ACMDistribution x))))

(extend-protocol density-function
  drec
  (pdf
    ([d x] (.density ACMDistribution x))
    ([d] (fn [x] (pdf ACMDistribution x))))
  (log-pdf
    ([d x] (.logDensity ACMDistribution x))
    ([d] (fn [x] (log-pdf ACMDistribution x)))))

(extend-protocol distribution-function
  drec
  (cdf
    ([d] (fn [x] (.cumulativeProbability ACMDistribution x)))
    ([d x] (.cumulativeProbability ACMDistribution x))))

(extend-protocol support
  drec
  (support-lower [d] (.getSupportLowerBound ACMDistribution))
  (support-upper [d] (.getSupportUpperBound ACMDistribution)))

(extend-protocol inverse-distribution-function
  drec
  (icdf
    ([d] (fn [x] (.inverseCumulativeProbability ACMDistribution x)))
    ([d x] (.inverseCumulativeProbability ACMDistribution x))))

(extend-protocol first-moment
  drec
  (mean [d] (.getNumericalMean ACMDistribution)))

(extend-protocol second-central-moment
  drec
  (variance [d] (.getNumericalVariance ACMDistribution)))

(extend-protocol random
  drec
  (sample
    ([d] (.sample ACMDistribution))
    ([d n] (take n (repeatedly #(sample ACMDistribution))))))
