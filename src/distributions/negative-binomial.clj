(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.PascalDistribution)

(defrecord NegativeBinomial [failures probability]
  random
  (sample [d] (.sample (new PascalDistribution failures probability)))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d x] (.probability (new PascalDistribution failures probability) x))
  (P [d] (fn [x] (P d x)))
  distribution-function
  (cdf [d x] (.cumulativeProbability (new PascalDistribution failures probability) x))
  (cdf [d] (fn [x] (cdf d x)))
  inverse-distribution-function
  (icdf [d x] (.inverseCumulativeProbability (new PascalDistribution failures probability) x))
  (icdf [d] (fn [x] (icdf d x)))
  mass-function
  (pmf [d] (fn [x] (pdf d x)))
  (pmf [d x] (if (integer? x) (.probability (new PascalDistribution failures probability) x) 0))
  (log-pmf [d x] (if (integer? x) (.logProbability (new PascalDistribution failures probability) x) Double/NEGATIVE_INFINITY))
  (log-pmf [d] (fn [x] (log-pmf d x)))
  density-function
  (pdf [d] (fn [x] (pdf d x)))
  (pdf [d x] (if (integer? x) (.probability (new PascalDistribution failures probability) x) 0))
  (log-pdf [d x] (if (integer? x) (.logProbability (new PascalDistribution failures probability) x) Double/NEGATIVE_INFINITY))
  (log-pdf [d] (fn [x] (log-pmf d x)))
  support
  (support-lower [d] (.getSupportLowerBound (new PascalDistribution failures probability)))
  (support-upper [d] (.getSupportUpperBound (new PascalDistribution failures probability)))
  first-moment
  (mean [d] (.getNumericalMean (new PascalDistribution failures probability)))
  second-central-moment
  (variance [d] (.getNumericalVariance (new PascalDistribution failures probability))))

(defn negative-binomial
  [failures probability]
  (new NegativeBinomial failures probability))

