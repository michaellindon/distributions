(in-ns 'distributions.core)

(defrecord Mixture [components probabilities]
  density-function
  (pdf [d x] (reduce + 0 (map (fn [c p] (* p (pdf c x))) components probabilities)))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d x]
    (let [logprob (map log probabilities)
          logcomp (map #(log-pdf % x) components)]
      (log-sum-exp (map + logprob logcomp))))
  (log-pdf [d] (fn [x] (log-pdf d x)))
  )


(defn mixture [components probabilities & {:keys [log?] :or {log? false}}]
  (let [prob (if log? (normalize-log probabilities) (normalize probabilities))]
    (Mixture. components prob))  )

(extend-protocol distribution-function
  Mixture
  (cdf
    ([d] (fn [x] (cdf d x)))
    ([d x] (reduce + 0 (map (fn [c p] (* p (cdf c x))) (:components d) (:probabilities d))))))

(extend-protocol support
  Mixture
  (support-lower [d] (reduce min (map support-lower (:components d))))
  (support-upper [d] (reduce max (map support-upper (:components d)))))

(extend-protocol first-moment
  Mixture
  (mean [d] (reduce + (map (fn [c p] (* p (mean c))) (:components d) (:probabilities d)))))

(extend-protocol second-central-moment
  Mixture
  (variance [d]
    (let [mu (mean d)]
      (map (fn [c p] (* p (+ (variance c) (square (- (mean c) mu))))) (:components d) (:probabilities d)))))

(extend-protocol random
  Mixture
  (sample
    ([d]
     (let [weights (:probabilities d)
           components (:components d)
           n (count weights)
           i (sample (discrete-integer (range 0 n) weights))]
       (sample (nth components i))))
    ([d n] (take n (repeatedly #(sample d))))))


(defmethod marginal [java.lang.Object distributions.core.Mixture]
  [likelihood {components :components probabilities :probabilities}]
  (mixture (map #(marginal likelihood %) components) probabilities))
