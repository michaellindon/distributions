(in-ns 'distributions.core)

(import org.apache.commons.math3.distribution.GammaDistribution)
(import org.apache.commons.math3.distribution.NormalDistribution)

(defrecord NormalInverseGamma [location lambda shape scale]
  random
  (sample [d] (let [s2 (sample (inverse-gamma shape scale))
                    mu (sample (normal location (/ s2 lambda)))]
                [mu s2]))
  (sample [d n] (take n (repeatedly #(sample d))))
  probability-function
  (P [d [mu s2]] (if (< s2 0)
                   0
                   (* (pdf (inverse-gamma shape scale) s2)
                      (pdf (normal location (/ s2 lambda)) mu))))
  (P [d] (fn [x] (pdf d x)))
  density-function
  (pdf [d [mu s2]] (if (< s2 0)
                     0
                     (* (pdf (inverse-gamma shape scale) s2)
                        (pdf (normal location (/ s2 lambda)) mu))))
  (pdf [d] (fn [x] (pdf d x)))
  (log-pdf [d [mu s2]] (if (< s2 0)
                         Double/NEGATIVE_INFINITY
                         (+ (log-pdf (inverse-gamma shape scale) s2)
                            (log-pdf (normal location (/ s2 lambda)) mu))))
  (log-pdf [d] (fn [x] (pdf d x))))

(defn normal-inverse-gamma [location prec shape scale] (new NormalInverseGamma location prec shape scale))

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
        post-shape (+ prior-shape (/ n 2))]
    (normal-inverse-gamma post-mean post-lambda post-shape post-scale)))

(defmethod marginal [distributions.core.Normal distributions.core.NormalInverseGamma]
  [likelihood prior]
  (let [{prior-shape :shape
         prior-lambda :lambda
         prior-scale :scale
         prior-location :location} prior]
    (t-distribution (* 2 prior-shape)
                    prior-location
                    (sqrt (/ (* prior-scale (+ 1 prior-lambda))
                             (* prior-lambda prior-shape))))))

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
  (def prior-shape 10)
  (def prior-scale 3)
  (- (log-pdf (posterior Y (normal :mu :s2) (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale)) [mu1 s21])
     (log-pdf (posterior Y (normal :mu :s2) (normal-inverse-gamma  prior-mean prior-lambda prior-shape prior-scale)) [mu2 s22]))

  (- (+ (log-pdf (normal prior-mean (/ s21 prior-lambda)) mu1)
        (log-pdf (inverse-gamma prior-shape prior-scale) s21)
        (reduce + (map #(log-pdf (normal mu1 s21) %) Y)))
     (+ (log-pdf (normal prior-mean (/ s22 prior-lambda)) mu2)
        (log-pdf (inverse-gamma prior-shape prior-scale) s22)
        (reduce + (map #(log-pdf (normal mu2 s22) %) Y))))

  (double (+ prior-scale (/ (+ (* (/ (* n prior-lambda) (+ n prior-lambda)) (square (- (mean Y) prior-mean)))
                               (reduce + (map square (sub Y (mean Y))))) 2)))

  (double (:scale (posterior Y (normal :mu :s2) (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale))))

  (variance (t-distribution 4))
  (variance (t-distribution 4 0 10))

  (* (variance (t-distribution (* 2 prior-shape)))  (/ (* prior-scale (+ 1 prior-lambda)) (* prior-lambda prior-shape)))
  (variance (map (fn [[mu s2]] (sample (normal mu s2))) (sample (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale) 100000)))
  (variance (marginal (normal :mu :s2) (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale)))
  (variance (sample (marginal (normal :mu :s2) (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale)) 10000))

  (* (variance (t-distribution (* 2 prior-shape))) (/ (* prior-lambda prior-shape) prior-scale))
  (variance (map first
                 (sample (normal-inverse-gamma prior-mean prior-lambda prior-shape prior-scale) 10000)))

  (double (variance (inverse-gamma 10 10)))
  (variance (sample (inverse-gamma 10 10) 100000))

  (def d (inverse-gamma 10 3))
  (double (variance d))
  (variance (sample d 100000))

  (variance (map (fn [[mu s2]] (sample (normal mu s2))) (sample (normal-inverse-gamma 4 2 1000 4) 100000)))
  (variance (map second (sample (normal-inverse-gamma 4 3 5 4) 10000)))
  (variance (inverse-gamma 5 4))
  (variance (sample (marginal (normal :mu :s2) (normal-inverse-gamma 4 2 1000 4)) 100000))
  (variance (marginal (normal :mu :s2) (normal-inverse-gamma 4 2 1000 4)))
  (variance (t-distribution 30000 3 10)))
