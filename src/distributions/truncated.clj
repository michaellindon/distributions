(in-ns 'distributions.core)

(defrecord TruncatedDistribution [distribution lower upper F-lower F-upper])
(defn truncated [distribution lower upper]
  (TruncatedDistribution. distribution lower upper (cdf distribution lower) (cdf distribution upper)))

(extend-protocol density-function
  TruncatedDistribution
  (pdf
    ([d x]
     (if (or (> x (:upper d)) (<= x (:lower d)))
       0
       (/ (pdf (:distribution d) x) (- (:F-upper d) (:F-lower d)))))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x]
     (if (or (> x (:upper d)) (<= x (:lower d)))
       Double/NEGATIVE_INFINITY
       (- (log-pdf (:distribution d) x) (log (- (:F-upper d) (:F-lower d))))))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol distribution-function
  TruncatedDistribution
  (cdf
    ([d] (fn [x] (.cumulativeProbability d x)))
    ([d x]
     (cond
       (> x (:upper d)) 1.0
       (<= x (:lower d)) 0.0
       :else (/ (- (cdf (:distribution d) x) (:F-lower d)) (- (:F-upper d) (:F-lower d)))))))

(extend-protocol support
  TruncatedDistribution
  (support-lower [d] (:lower d))
  (support-upper [d] (:upper d)))

(extend-protocol inverse-distribution-function
  TruncatedDistribution
  (icdf
    ([d] (fn [x] (icdf d x)))
    ([d x]
     (let [base (:distribution d)
           Fa (:F-lower d)
           Fb (:F-upper d)]
       (icdf base (+ (* x (- Fb Fa)) Fa))))))
