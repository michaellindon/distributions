(in-ns 'distributions.core)

(defrecord InverseGaussianDistribution [mean shape])
(defn inverse-gaussian
  [mean shape]
  (InverseGaussianDistribution. mean shape))

(extend-protocol random
  InverseGaussianDistribution
  (sample
    ([d]
     (let [shape (:shape d)
           mean (:mean d)
           v (sample (normal 0 1))
           y (square v)
           x (+ mean
                (* 0.5 (square mean) (/ y shape))
                (negate
                 (* 0.5
                    (/ mean shape)
                    (sqrt (+ (* 4 mean shape y) (* (square mean) (square y)))))))
           z (sample (uniform 0 1))]
       (if (< z (/ mean (+ mean x))) x (/ (square mean) x))))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol distribution-function
  InverseGaussianDistribution
  (cdf
    ([d] (fn [x] (cdf d x)))
    ([d x]
     (if (<= x 0)
       0
       (let [shape (:shape d)
             mean (:mean d)]
         (+
          (probit (* (sqrt (/ shape x)) (dec (/ x mean))))
          (*
           (exp (* 2 (/ shape mean)))
           (probit (negate (* (sqrt (/ shape x)) (inc (/ x mean))))))))))))

(extend-protocol density-function
  InverseGaussianDistribution
  (pdf
    ([d x] (exp (log-pdf d x)))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x]
     (if (<= x 0)
       Double/NEGATIVE_INFINITY
       (let [shape (:shape d)
             mean (:mean d)]
         (* 0.5 (+ (log shape)
                   (negate (log (* 2 Math/PI x x x)))
                   (negate
                    (/
                     (* shape (square (- x mean)))
                     (* mean mean x))))))))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol probability-function
  InverseGaussianDistribution
  (P
    ([d] (fn [x] 0))
    ([d x] 0)
    ))

(extend-protocol support
  InverseGaussianDistribution
  (support-lower [d] 0)
  (support-upper [d] Double/POSITIVE_INFINITY)
  )

(extend-protocol first-moment
  InverseGaussianDistribution
  (mean [d] (:mean d)))

(extend-protocol second-central-moment
  InverseGaussianDistribution
  (variance [d] (/ (pow (:mean d) 3) (:shape d))))
