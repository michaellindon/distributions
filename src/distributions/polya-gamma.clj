(in-ns 'distributions.core)

(defn jacobi-a
  ([x n z] (* (cosh z) (exp (negate (* z z 0.5 x))) (jacobi-a x n)))
  ([x n]
   (if (< x 0.64)
     (* Math/PI (+ n 0.5) (pow (/ 2 (* x Math/PI)) 1.5) (exp (negate (/ (* 2 (square (+ n 0.5))) x))))
     (* Math/PI (+ n 0.5) (exp (negate (* 0.5 x (square (* Math/PI (+ n 0.5)))))) )))
  )

(defn jacobi-iterator [x z [sp ip]]
  (let [i (inc ip)
        s (if (even? i) (+ sp (jacobi-a x i z)) (- sp (jacobi-a x i z)) )]
    [s i]))


(defn randJ [z]
  (let [t 0.64
        ex (exponential (+ (* 0.5 z z) (/ (square Math/PI) 8)))
        ig (inverse-gaussian (inv (abs z)) 1)
        p (* (inc (exp (negate (* 2 (abs z))))) (cdf ig t))
        q (* (cosh z) 0.5 Math/PI (mean ex) (- 1 (cdf ex t)))
        ig-prob (/ p (+ p q))]
    (loop []
      (let [g (mixture [(truncated ig 0 t) (truncated ex t Double/POSITIVE_INFINITY)] [ig-prob (- 1 ig-prob)])
            X (sample g)
            U (sample (uniform 0 (* (+ p q) (pdf g X))))
            stream (drop 1 (iterate (partial jacobi-iterator X z) [(jacobi-a X 0 z) 0]))
            S (first (drop-while (fn [[sn n]] (not (if (odd? n) (< U sn) (> U sn)))) stream))]
        (if (odd? (second S)) X (recur))
        )
      )
    ))

(defn not-converged? [[[sn _] [sp _]]] (> (abs (- sn sp)) (java.lang.Math/ulp 1.0)))

(defn jacobi-density [z x]
  (let [stream (iterate (partial jacobi-iterator x z) [(jacobi-a x 0 z) 0])
        zip-stream (map vector (drop 1 stream) stream)]
    (first (first (first (drop-while not-converged? zip-stream))))))

(defrecord PolyaGammaDistribution [z])
(defn polya-gamma
  [z]
  (PolyaGammaDistribution. z))

(extend-protocol probability-function
  PolyaGammaDistribution
  (P
    ([d] (fn [x] 0))
    ([d x] 0)
    ))

(extend-protocol density-function
  PolyaGammaDistribution
  (pdf
    ([d x] (* 4 (jacobi-density (* 0.5 (:z d)) (* 4 x))))
    ([d] (fn [x] (pdf d x))))
  (log-pdf
    ([d x] (log (pdf d x)))
    ([d] (fn [x] (log-pdf d x)))))

(extend-protocol random
  PolyaGammaDistribution
  (sample
    ([d] (* 0.25 (randJ (* 0.5 (:z d)))))
    ([d n] (take n (repeatedly #(sample d))))))



(comment
(defn naiverand [z]
  (/ (reduce + (for [i (range 1 100)]
                 (/ (sample (gamma 1 1)) (+ (square (/ z (* 2 Math/PI))) (square (- i 0.5)))))) (* 2 Math/PI Math/PI))))
