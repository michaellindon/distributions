(in-ns 'distributions.core)

(defrecord MixtureDistribution [components probabilities])
(defn mixture [components probabilities]
  (MixtureDistribution. components probabilities))

(extend-protocol density-function
  MixtureDistribution
  (pdf
    ([d x] (reduce + 0 (map (fn [c p] (* p (pdf c x))) (:components d) (:probabilities d))))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (log (pdf d x)))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol distribution-function
  MixtureDistribution
  (cdf
    ([d] (fn [x] (.cumulativeProbability d x)))
    ([d x] (reduce + 0 (map (fn [c p] (* p (cdf c x))) (:components d) (:probabilities d))))))

(extend-protocol support
  MixtureDistribution
  (support-lower [d] (reduce min (map support-lower (:components d))))
  (support-upper [d] (reduce max (map support-upper (:components d)))))

(extend-protocol first-moment
  MixtureDistribution
  (mean [d] (reduce + (map (fn [c p] (* p (mean c))) (:components d) (:probabilities d)))))

(extend-protocol second-central-moment
  MixtureDistribution
  (variance [d]
    (let [mu (mean d)]
      (map (fn [c p] (* p (+ (variance c) (square (- (mean c) mu))))) (:components d) (:probabilities d)))))

(extend-protocol random
  MixtureDistribution
  (sample
    ([d]
     (let [weights (:probabilities d)
           components (:components d)
           n (count weights)
           i (sample (discrete-integer (range 0 n) weights))]
       (sample (get components i))))
    ([d n] (take n (repeatedly #(sample d))))))
