(in-ns 'distributions.core)


(defrecord MultivariateNormalDistribution [mean variance])
(defn mvnormal
  ([mean variance]
   (MultivariateNormalDistribution. mean variance))
  ([variance]
   (let [n (column-count variance)
         mean (zero-vector n)]
     (mvnormal mean variance))))
(extend-protocol density-function
  MultivariateNormalDistribution
  (log-pdf
    ([d x]
     (let [mean-vector (:mean d)
           cov-matrix (:variance d)
           n (ecount x)
           residual (sub x mean-vector)
           exponent (* 0.5 (dot residual (mmul (inverse cov-matrix) residual)))
           normalizer (+ (* 0.5 (log (det cov-matrix)))
                         (* n 0.5 (log (* 2 Math/PI))))]
       (negate (+ normalizer exponent))
       ))
    ([d] (fn [x] (log-pdf d x))))
  (pdf
    ([d x] (exp (log-pdf d x)))
    ([d] (fn [x] (pdf d x)))))
(extend-protocol first-moment
  MultivariateNormalDistribution
  (mean [d] (:mean d)))
(extend-protocol random
  MultivariateNormalDistribution
  (sample
    ([d]
     (let [cov-matrix (:variance d)
           mean-vector (:mean d)
           ncols (column-count cov-matrix)
           z (sample (normal 0 1) ncols)
           d (la/svd cov-matrix)
           {U :U S :S V* :V*} (la/svd cov-matrix)
           D (diagonal-matrix (map (fn [x] (sqrt (max x 0))) S))]
       (add mean-vector (mmul U D z))))
    ([d n] (take n (repeatedly #(sample d))))))
