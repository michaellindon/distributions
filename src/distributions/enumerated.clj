(in-ns 'distributions.core)


(defrecord Enumerated [items probabilities discrete imap])

(defn enumerated [items probabilities]
  (let [integers (range 0 (count items))
        Z (reduce + probabilities)
        prob (map (fn [x] (/ x Z)) probabilities)]
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
