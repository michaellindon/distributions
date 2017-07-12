(in-ns 'distributions.core)

(defn mh-step
  ([target proposal x-old] (mh-step target proposal x-old false))
  ([target proposal x-old log-target?]
   (let [x-new (sample (proposal x-old))
         log-target (if log-target? target #(log (target %)))
         acc-prob (min 1 (exp (+ (log-target x-new)
                                 (negate (log-target x-old))
                                 (log-pdf (proposal x-new) x-old)
                                 (negate (log-pdf (proposal x-old) x-new)))))]
     (if (= 1 (sample (bernoulli acc-prob))) x-new x-old))))

(defn metropolis-hastings
  ([target proposal x-old] (metropolis-hastings target proposal x-old false))
  ([target proposal x-old log-target?]
   (let [mh-iterator #(mh-step target proposal % log-target?)]
     (drop 1 (iterate mh-iterator x-old)))))

(defn slice-step
  "Sample a univariate unnormalized log-density g from position x with length L"
  [g w x]
  (let [y (+ (g x) (negate (sample (exponential 1))))
        u (sample (uniform 0 1))
        lower-bound (first (filter (fn [x] (< (g x) y))
                                   (iterate  (fn [x] (- x w))
                                             (- x (* u w)))))
        upper-bound (first (filter (fn [x] (< (g x) y))
                                   (iterate (fn [x] (+ x w))
                                            (+ x (* (- 1 u) w)))))]
    (loop [l lower-bound
           u upper-bound]
      (let [z (sample (uniform l u))]
        (cond
          (> (g z) y) z
          (> z x) (recur l z)
          :else (recur z u))))))

(defn slice-sampler
  ([target width seed]
   (slice target width seed false))
  ([target width seed log-target?]
   (let [log-target (if log-target? target #(log (target %)))
         slice-iterator (partial slice-step log-target width)]
     (drop 1 (iterate slice-iterator seed)))))

(defn accept-reject [f g c]
  (let [x (repeatedly #(sample g))
        u (map #(sample (uniform 0 (* c (pdf g %)) )) x)
        ux (map vector x u)]
    (map first (filter (fn [[x u]] (< u (f x))) ux))))
