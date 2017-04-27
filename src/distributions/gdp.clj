(in-ns 'distributions.core)

(defrecord GeneralizedDoubleParetoDistribution [scale shape])
(defn gdp
  [scale shape]
  (GeneralizedDoubleParetoDistribution. scale shape))

(extend-protocol probability-function
  GeneralizedDoubleParetoDistribution
  (P
    ([d] (fn [x] 0))
    ([d x] 0)
    ))

(extend-protocol density-function
  GeneralizedDoubleParetoDistribution
  (pdf
    ([d x] (exp (log-pdf d x)))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x]
     (let [shape (:shape d)
           scale (:scale d)]
       (- (log 0.5) (log scale) (* (inc shape) (log (inc (/ (abs x) (* shape scale)) ))))))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol proximal
  GeneralizedDoubleParetoDistribution
  (prox
    ([d]
     (let [shape (:shape d)
           scale (:scale d)]
       (fn [h x]
         (let [g (* h (inc shape))
               a (* shape scale)
               d (max 0 (- (* a (abs x)) g))]
           (* 0.5 (signum x) (+ (abs x) (negate a) (sqrt (+ (* 4 d) (square (- a (abs x)))))))))))
    ([d h x]
     ((prox d) h x))))

(extend-protocol random
  GeneralizedDoubleParetoDistribution
  (sample
    ([d]
     (let [shape (:shape d)
           scale (:scale d)
           eta (* shape scale)
           lambda (sample (gamma shape eta))
           tau (sample (exponential (/ (square lambda) 2)))]
       (sample (normal 0 tau))))
    ([d n] (take n (repeatedly #(sample d))))))
