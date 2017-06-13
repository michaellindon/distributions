(in-ns 'distributions.core)

(defrecord DirichletDistribution [alpha])

(defn dirichlet [alpha] (DirichletDistribution. alpha))

(extend-protocol random
  DirichletDistribution
  (sample
    ([d] (let [components (map #(sample (gamma % 1)) (:alpha d))
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
     (let [alphas (:alpha d)
           a-sum (reduce + alphas)
           log-gammas (map  log-gamma-fn alphas)
           log-B (- (reduce + log-gammas) (log-gamma-fn a-sum))
           ]
       (- (reduce + (map (fn [y a] (* (dec a) (Math/log y))) x alphas)) log-B)))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol first-moment
  DirichletDistribution
  (mean
    [d]
    (let [alphas (:alpha d)
          a-sum (reduce + alphas)]
      (mapv #(/ % a-sum) alphas))))

(extend-protocol second-central-moment
  DirichletDistribution
  (variance [d]
    (let [alphas (:alpha d)
          a0 (reduce + alphas)
          denom (* (inc a0) (square a0))
          d (count alphas)]
      (reshape (for [i (range 0 d)
                     j (range 0 d)]
                 (if (= i j)
                   (/ (* (nth alphas i) (- a0 (nth alphas i))))
                   (/ (negate (* (nth alphas i) (nth alphas j))) denom))) [d d]))))
