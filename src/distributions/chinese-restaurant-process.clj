(in-ns 'distributions.core)


(defrecord ChineseRestaurantProcess [colours concentration discount])

(defn cr-process
  ([concentration discount] (cr-process {} concentration discount))
  ([colours concentration discount] (ChineseRestaurantProcess. colours concentration discount)))

(defrecord CRPRealization [colours concentration discount])

(defn crp-realization [colours concentration discount] (CRPRealization. colours concentration discount))


(extend-protocol random
  ChineseRestaurantProcess
  (sample
    ([d] (let [colours (atom (:colours d))
               {concentration :concentration
                discount :discount} d]
           (crp-realization colours concentration discount)))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol random
  CRPRealization
  (sample
    ([d] (let [{concentration :concentration
                discount :discount
                colour-ref :colours} d
               colours @colour-ref
               K (count colours)
               n (reduce + (vals colours))]
           (if (empty? colours)
             (do
               (swap! colour-ref assoc 1 1)
               1)
             (let [repeat-prob (map (fn [[c nc]] (/ (- nc discount) (+ concentration n) )) colours)
                   new-prob (/ (+ concentration (* discount K)) (+ n concentration))
                   proposal (discrete-integer (conj (keys colours) (inc K)) (conj repeat-prob new-prob))
                   newc (sample proposal)]
               (do
                 (if (contains? colours newc)
                   (swap! colour-ref assoc newc (inc (get colours newc)))
                   (swap! colour-ref assoc newc 1))
                 newc))
             )))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol first-moment
  CRPRealization
  (mean
    [d]
    (let [{concentration :concentration
           discount :discount} d]
      (fn [n] (- (exp (+
                       (log-gamma-fn (+ concentration discount n))
                       (log-gamma-fn (inc concentration))
                       (negate (+
                                (log discount)
                                (log-gamma-fn (+ concentration n))
                                (log-gamma-fn (+ concentration discount) )))))
                 (/ concentration discount))))))

(comment
  ((mean (sample (cr-process 4.3 0.5))) 100)
  (def foo (sample (cr-process 4.3 0.5) 10000))
  (def bar (map #(sample % 100) foo))
  (double (mean (map (comp count distinct) bar))))
