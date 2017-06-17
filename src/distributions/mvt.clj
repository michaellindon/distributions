(in-ns 'distributions.core)


(defrecord MultivariateT [df mean covariance])
(defn mvt
  ([df mean covariance]
   (MultivariateT. df mean covariance))
  ([df covariance]
   (let [n (column-count covariance)
         mean (zero-vector n)]
     (mvt mean covariance))))
(extend-protocol density-function
  MultivariateT
  (log-pdf
    ([d x]
     (let [df (:df d)
           mean-vector (:mean d)
           cov-matrix (:covariance d)
           n (ecount x)
           residual (sub x mean-vector)
           kernel (* (negate (* 0.5 (+ df n))) (log (inc (/ (dot residual (mmul (inverse cov-matrix) residual)) df))))
           normalizer (- (log-gamma-fn (* 0.5 (+ df n))) (log-gamma-fn (* 0.5 df)) (* 0.5 n (log (* Math/PI df))) (* 0.5 (log (det cov-matrix))))]
       (+ normalizer kernel)))
    ([d] (fn [x] (log-pdf d x))))
  (pdf
    ([d x] (exp (log-pdf d x)))
    ([d] (fn [x] (pdf d x)))))
(extend-protocol first-moment
  MultivariateT
  (mean [d] (:mean d)))
(extend-protocol random
  MultivariateT
  (sample
    ([d]
     (let [df (:df d)
           cov-matrix (:covariance d)
           mean-vector (:mean d)
           Z (sample (mvnormal cov-matrix))
           u (sample (gamma (* 0.5 df) (* 0.5 df)))]
       (add mean-vector (div Z (sqrt u)))))
    ([d n] (take n (repeatedly #(sample d))))))
