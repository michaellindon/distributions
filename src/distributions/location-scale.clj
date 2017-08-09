(in-ns 'distributions.core)

(defrecord LocationScaleDistribution [distribution location scale]
  support
  (support-lower [d] (+ location (support-lower distribution)))
  (support-upper [d] (+ location (support-upper distribution)))
  first-moment
  (mean [d] (+ location (mean distribution)))
  second-central-moment
  (variance [d] (* scale scale (variance distribution))))

(defn location-scale
  [distribution location scale]
  (LocationScaleDistribution. distribution location scale))

(extend-protocol density-function
  LocationScaleDistribution
  (log-pdf
    ([d] (fn [x] (log-pdf d x)))
    ([d x] (- (log-pdf (:distribution d) (/ (- x (:location d)) (:scale d))) (log (:scale d)))))
  (pdf
    ([d] (fn [x] (pdf d x)))
    ([d x] (/ (pdf (:distribution d) (/ (- x (:location d)) (:scale d))) (:scale d)))))

(extend-protocol distribution-function
  LocationScaleDistribution
  (cdf
    ([d] (fn [x] (cdf d x)))
    ([d x] (cdf (:distribution d) (/ (- x (:location d)) (:scale d))))))

(extend-protocol random
  LocationScaleDistribution
  (sample
    ([d] (+ (:location d) (* (:scale d) (sample (:distribution d)))))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol inverse-distribution-function
  LocationScaleDistribution
  (icdf
    ([d] (fn [x] (icdf d x)))
    ([d x] (+ (:location d) (* (:scale d) (icdf (:distribution d) x))))))

