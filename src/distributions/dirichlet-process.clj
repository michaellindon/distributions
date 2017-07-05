(ns distributions.dirichlet-process
  (:require [clojure.set :as set]))
(use 'distributions.core)


(defrecord DirichletProcess [concentration base-measure])

(defn dirichlet-process [concentration base-measure] (DirichletProcess. concentration base-measure))

(defrecord DPRealization [points-ref concentration base-measure])

(defn dp-realization [points-ref concentration base-measure] (DPRealization. points-ref concentration base-measure))

(extend-protocol random
  DirichletProcess
  (sample
    ([d] (let [points (atom '())
               {concentration :concentration
                base-measure :base-measure} d]
           (dp-realization points concentration base-measure)))
    ([d n] (take n (repeatedly #(sample d))))))

(extend-protocol random
  DPRealization
  (sample
    ([d] (let [{points-ref :points-ref
                concentration :concentration
                base-measure :base-measure} d
               points @points-ref
               n (count points)]
           (if (zero? n)
             (let [newval (sample base-measure)]
               (do
                 (swap! points-ref conj newval)
                 newval))
             (let [disc-prob (/ 1 (+ n concentration))
                   base-prob (/ concentration (+ n concentration))
                   disc-measure (discrete-real points (take n (repeat (/ 1 n))))
                   mix (mixture [disc-measure base-measure] [disc-prob base-prob])
                   newval (sample mix)]
               (do
                 (swap! points-ref conj newval)
                 newval))
             )))
    ([d n] (take n (repeatedly #(sample d))))))


(ns distributions.core)
(def concentration 0.0001)
(def G (normal 0 3))
(defn posterior [y G]
  (let [n (count y)
        s2 (variance G)
        ybar (mean y)
        postvar (/ 1 (+ (/ 1 s2) n))]
    (normal (* postvar ybar) postvar)))

(defn marginal [F G]
  (normal (mean G) (+ (variance F) (variance G))))

(defn positions
  [pred coll]
  (keep-indexed (fn [idx x]
                  (when (pred x)
                    idx))
                coll))

(def y (sample (mixture [(normal 0 1) (normal -1 1) (normal 1 1)] [(/ 1 3) (/ 1 3) (/ 1 3)]) 100))
(defn update-c [observations likelihood base-measure labels]
  (let [n (count labels)
        prior-predictive (marginal likelihood base-measure)]
    (loop [cs labels
           i 0]
      (if (== i n)
        cs
        (let [freqs (dissoc (frequencies cs) (get cs i))
              c-vals (keys freqs)
              c-counts (vals freqs)
              y (get observations i)
              c-probs (map (fn [[c nc]]
                             (let [idx (positions #(identical? c %) cs)
                                   y-ic (for [j idx] (get observations j))
                                   posterior-ic (posterior y-ic base-measure)
                                   posterior-predictive (marginal likelihood posterior-ic)]
                               (* (/ nc (+ n -1 concentration)) (pdf posterior-predictive y)))) freqs)
              not-c-prob (* (/ concentration (+ n -1 concentration)) (pdf prior-predictive y))
              proposal (enumerated (conj c-vals (gensym)) (conj c-probs not-c-prob))]
          (recur (assoc cs i (sample proposal)) (inc i))
          ))
      )))

(def foo (doall (sample (normal 0 1) 100)))
(def draws (take 100 (iterate (partial update-c (into [] y) (normal 0 1) (normal 0 3)) (into [] foo))))
(map (comp count distinct) draws)
