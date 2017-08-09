(in-ns 'distributions.core)


(defrecord Enumerated [items probabilities discrete imap])

(defn enumerated [items probabilities & {:keys [log?]
                                         :or {log? false}}]
  (let [integers (range 0 (count items))
        prob (if log? (normalize-log probabilities) (normalize probabilities))]
    (Enumerated. items prob (discrete-integer integers prob) (into {} (map vector items integers )))))

(extend-protocol random
  Enumerated
  (sample
    ([d] (let [{discrete :discrete
                items :items} d]
           (nth items (sample discrete))))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol probability-function
  Enumerated
  (P
    ([d x]
     (let [{imap :imap
            discrete :discrete} d]
       (pmf discrete (imap x))))
    ([d] (fn [x] (P d x)))))
