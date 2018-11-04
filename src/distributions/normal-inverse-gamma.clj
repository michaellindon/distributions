(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.GammaDistribution)
(import org.apache.commons.math3.distribution.NormalDistribution)

(defrecord NormalInverseGamma [location lambda shape scale]
  random
  (sample [d] (let [s2 (/ 1 (.sample (new GammaDistribution shape (/ 1 scale))))
                    mu (.sample (new NormalDistribution location (sqrt (/ s2 lambda))))]
                [mu s2]))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d [mu s2]] (* (pdf (inverse-gamma shape scale) s2)
                      (pdf (normal location (/ s2 lambda)) mu)))
  (P [d] (fn [x] (pdf d x)))
  density-function
  (pdf [d [mu s2]] (* (pdf (inverse-gamma shape scale) s2)
                      (pdf (normal location (/ s2 lambda)) mu)))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d [mu s2]] (+ (log-pdf (inverse-gamma shape scale) s2)
                          (log-pdf (normal location (/ s2 lambda)) mu)))
  (log-pdf [d] (fn [x] (pdf d x))))

(defmethod posterior [distributions.core.Normal distributions.core.NormalInverseGamma]
  [y likelihood prior]
  (let [prior-mean (:location prior)
        prior-lambda (:lambda prior)
        prior-shape (:shape prior)
        prior-scale (:scale prior)
        ybar (mean y)
        n (count y)
        post-lambda (+ n prior-lambda)
        post-mean (/ (+ (* n ybar)
                        (* prior-lambda prior-mean))
                     post-lambda)
        post-scale (+ prior-scale (/ (+
                                      (reduce + (map square (sub y post-mean)))
                                      (* prior-lambda
                                         (square (- post-mean
                                                    prior-mean))))
                                     2))
        post-shape (+ prior-shape (/ n 2))
        ]
    (normal-inverse-gamma post-mean post-lambda post-shape post-scale)))

(defn normal-inverse-gamma [location prec shape scale] (new NormalInverseGamma location prec shape scale))

(defmethod marginal [distributions.core.Normal distributions.core.NormalInverseGamma]
  [likelihood prior]
  (let [{prior-shape :shape
         prior-lambda :lambda
         prior-scale :scale
         prior-location :location} prior]
    (t-distribution (* 2 prior-shape)
                    prior-location
                    (/ (* prior-scale (+ 1 prior-lambda))
                       (* prior-lambda prior-shape)))))

(comment 
(sample (marginal (normal :mu :s2) (normal-inverse-gamma 2 3 1 2)))

(def Y [1 1 1 2 2 3 1 2 5])
(def n (count Y))

(def mu1 3)
(def s21 10)
(def mu2 3)
(def s22 1)
(def prior-mean 3)
(def prior-lambda 2)
(def prior-shape 5)
(def prior-scale 3)
(- (log-pdf (posterior Y (normal :mu :s2) (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale)) [mu1 s21])
   (log-pdf (posterior Y (normal :mu :s2) (normal-inverse-gamma  prior-mean prior-lambda prior-shape prior-scale)) [mu2 s22]) )

(- (+ (log-pdf (normal prior-mean (/ s21 prior-lambda)) mu1)
      (log-pdf (inverse-gamma prior-shape prior-scale) s21)
      (reduce + (map #(log-pdf (normal mu1 s21) %) Y)))
   (+ (log-pdf (normal prior-mean (/ s22 prior-lambda)) mu2)
      (log-pdf (inverse-gamma prior-shape prior-scale) s22)
      (reduce + (map #(log-pdf (normal mu2 s22) %) Y)))
   )


(double (+ prior-scale (/ (+ (* (/ (* n prior-lambda) (+ n prior-lambda)) (square (- (mean Y) prior-mean)))
                            (reduce + (map square (sub Y (mean Y))))) 2)))

(double (:scale (posterior Y (normal :mu :s2) (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale))))
)
