(in-ns 'distributions.core)

(defrecord DirichletDistribution [concentration])

(defn dirichlet [concentration] (DirichletDistribution. concentration))

(extend-protocol random
  DirichletDistribution
  (sample
    ([d] (let [components (map #(sample (gamma % 1)) (:concentration d))
               component-sum (reduce + components)]
           (mapv #(/ % component-sum) components)))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol density-function
  DirichletDistribution
  (pdf
    ([d x] (Math/exp (log-pdf d x)))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x]
     (let [concentrations (:concentration d)
           a-sum (reduce + concentrations)
           log-gammas (map  log-gamma-fn concentrations)
           log-B (- (reduce + log-gammas) (log-gamma-fn a-sum))
           ]
       (- (reduce + (map (fn [y a] (* (dec a) (Math/log y))) x concentrations)) log-B)))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol first-moment
  DirichletDistribution
  (mean
    [d]
    (let [concentrations (:concentration d)
          a-sum (reduce + concentrations)]
      (mapv #(/ % a-sum) concentrations))))

(extend-protocol second-central-moment
  DirichletDistribution
  (variance [d]
    (let [concentrations (:concentration d)
          a0 (reduce + concentrations)
          denom (* (inc a0) (square a0))
          d (count concentrations)]
      (reshape (for [i (range 0 d)
                     j (range 0 d)]
                 (if (= i j)
                   (/ (* (nth concentrations i) (- a0 (nth concentrations i))))
                   (/ (negate (* (nth concentrations i) (nth concentrations j))) denom))) [d d]))))
