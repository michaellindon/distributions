(in-ns 'distributions.core)

(defrecord DirichletProcess [concentration base-measure])

(defn dirichlet-process [concentration base-measure] (DirichletProcess. concentration base-measure))

(defmethod posterior [clojure.lang.Keyword distributions.core.DirichletProcess]
  [data likelihood {concentration :concentration base-measure :base-measure}]
  (let [n (count data)
        disc-prob (/ 1 (+ n concentration))
        base-prob (/ concentration (+ n concentration))
        freqs (frequencies data)
        disc-measure (discrete-real (keys freqs) (vals freqs))]
    (dirichlet-process (+ n concentration) (mixture [disc-measure base-measure] [disc-prob base-prob]))))

(defmethod marginal [java.lang.Object distributions.core.DirichletProcess]
  [likelihood {concentration :concentration base-measure :base-measure}]
  base-measure)

(defn dp-update-labels [observations likelihood base-measure concentration labels]
  (let [n (count labels)
        possible-integers (set (range 0 n))
        prior-predictive (marginal likelihood base-measure)]
    (loop [cs labels
           i 0]
      (if (== i n)
        cs
        (let [c-i (remove-at cs i)
              yi (nth observations i)
              y-i (remove-at observations i)
              freqs (frequencies c-i)
              c-vals (keys freqs)
              alternative-integers (set/difference possible-integers (into #{} c-vals))
              c-probs (map (fn [[c nc]]
                             (let [idx (positions #(identical? c %) c-i) ;indices of c-i equal to c
                                   y-ic (for [j idx] (nth y-i j))
                                   posterior-predictive (posterior-predictive y-ic likelihood base-measure)]
                               (* (/ nc (+ (dec n) concentration)) (pdf posterior-predictive yi)))) freqs)
              not-c-prob (* (/ concentration (+ (dec n) concentration)) (pdf prior-predictive yi))
              proposal (enumerated (conj c-vals (first alternative-integers)) (conj c-probs not-c-prob))]
          (recur (assoc cs i (sample proposal)) (inc i))
          ))
      )))

(defn dp-update-params [observations likelihood base-measure concentration labels]
  (let [label-set (distinct labels)
        n (count observations)
        indices (range 0 n)
        index-sets (map (fn [x] (positions (fn [y] (identical? x y)) labels)) label-set)
        observation-sets (for [i index-sets] (for [j i] (nth observations j)))
        param-set (map (fn [obs] (sample (posterior obs likelihood base-measure))) observation-sets)
        label-param-map  (zipmap label-set param-set)]
    (for [c labels] (get label-param-map c))))


(defn rn-alg3 [observations likelihood base-measure concentration niter]
  (let [n (count observations)
        initial-labels (into [] (range 0 n))
        label-draws (take 100 (iterate (partial dp-update-labels observations likelihood base-measure concentration) initial-labels))
        param-draws (map (partial dp-update-params observations likelihood base-measure concentration ) label-draws)
        pred-fn-draws (map (fn [params] (pdf (marginal (normal :mu 1) (posterior-predictive params :G (dirichlet-process concentration base-measure))))) param-draws)]
    {:labels label-draws :parameters param-draws :densities pred-fn-draws }))


