(ns distributions.integration)

(defn expectation-mc [d n]
  (let [draws (doall (sample d n))]
    (fn [f]
      (->> draws (map f) mean))))

(defn expectation-qi [d n]
  (let [delta (/ 1.0 n)
        grid (doall (map (icdf d) (range delta 1 delta)))]
    (fn [f]
      (->> grid (map f) mean))))

(defn expectation
  ([d & {:keys [n method]
          :or {n 10000
               method "monte-carlo"}}]
   (case method
     "monte-carlo" (expectation-mc d n)
     "quantile-integration" (expectation-qi d n)
     "Invalid method")))

(defn kullback-leibler [f & {:keys [n method]
                             :or {n 10000
                                  method "monte-carlo"}}]
  (let [Ef (expectation f :n n :method method)
        Ef-log-f (Ef #(log-pdf f %))]
    (fn [g]
      (let [Ef-log-g (Ef #(log-pdf g %))]
        (- Ef-log-f Ef-log-g)))))
