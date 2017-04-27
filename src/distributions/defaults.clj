(in-ns 'distributions.core)

(extend-protocol random
  Object
  (sample
    ([d] (icdf d (sample (uniform 0 1))))
    ([d n] (take n (repeatedly #(sample d))))))
(extend-protocol inverse-distribution-function
  Object
  (icdf
    ([d] (fn [x] (icdf d x)))
    ([d x]
     (let [f #(- (cdf d %) x)
           w 1 ;(std d)
           m 0 ;(mean d)
           a (max (support-lower d) (first (drop-while #(> (f %) 0) (iterate #(- % w) m))))
           b (min (support-upper d) (first (drop-while #(< (f %) 0) (iterate #(+ % w) m))))]
       (bisection f a b)))))
(extend-protocol standard-deviation
  Object
  (std [d] (sqrt (variance d))))
(extend-protocol first-moment
  clojure.lang.Sequential
  (mean [coll] (if (empty? coll)
                 nil
                 (let [sum (reduce + coll)
                       n (count coll)]
                   (/ sum n)))))
(extend-protocol second-central-moment
  clojure.lang.Sequential
  (variance [coll] (if (empty? coll)
                     nil
                     (let [mu (mean coll)]
                       (mean (map (fn [x] (square (- x mu))) coll))))))
(extend-protocol random
  clojure.lang.LazySeq
  (sample
    ([d] (first (take 1 d)))
    ([d n] (take n d))))
