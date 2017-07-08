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

(defn positions
  [pred coll]
  (keep-indexed (fn [idx x]
                  (when (pred x)
                    idx))
                coll))

(defn remove-at [coll idx] (for [i (remove #(== % idx) (range 0 (count coll)))] (nth coll i)))


(defn update-labels [observations likelihood base-measure concentration labels]
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

(defn update-params [observations likelihood base-measure concentration labels]
  (let [label-set (distinct labels)
        n (count observations)
        indices (range 0 n)
        index-sets (map (fn [x] (positions (fn [y] (identical? x y)) labels)) label-set)
        observation-sets (for [i index-sets] (for [j i] (nth observations j)))
        param-set (map (fn [obs] (sample (posterior obs likelihood base-measure))) observation-sets)
        label-param-map  (zipmap label-set param-set)]
    (for [c labels] (get label-param-map c))))


(defmethod posterior [clojure.lang.Keyword distributions.dirichlet_process.DirichletProcess]
  [data likelihood {concentration :concentration base-measure :base-measure}]
  (let [n (count data)
        disc-prob (/ 1 (+ n concentration))
        base-prob (/ concentration (+ n concentration))
        freqs (frequencies data)
        disc-measure (discrete-real (keys freqs) (vals freqs))]
    (dirichlet-process (+ n concentration) (mixture [disc-measure base-measure] [disc-prob base-prob]))))

(defmethod marginal [java.lang.Object distributions.dirichlet_process.DirichletProcess]
  [likelihood {concentration :concentration base-measure :base-measure}]
  base-measure)



;Normal Mixture Model Example
(def likelihood (normal :mu 1))
(def base-measure (normal 0 4))
(def observations y)
(def concentration 0.1)
(def label-draws (take 100 (iterate (partial update-labels observations likelihood base-measure concentration) initial-labels)))
(def param-draws (map (partial update-params observations likelihood base-measure concentration ) label-draws))
(def pred-fn-draws (map (fn [params] (pdf (marginal (normal :mu 1) (posterior-predictive params :G (dirichlet-process concentration base-measure))))) param-draws))
(def seed (plot/function-plot (first pred-fn-draws) -5 5))
(def allplots (reduce (fn [chart newfun] (plot/add-function chart newfun -5 5)) seed (rest pred-fn-draws)))
(ic/view allplots)
(ic/view ((fn [chart newfun] (plot/add-function chart newfun -5 5)) seed (second pred-fn-draws)))

(def foo (plot/function-plot (pdf (marginal (normal :mu 1) (posterior-predictive (update-params y (normal :mu 1) (normal 0 4) 1 initial-labels) :G (dirichlet-process 1 (normal 0 4))))) -4 4))
(ic/view (plot/add-function foo (fn [x] (Math/sin x)) -5 5))
( (update-params y (normal :mu 1) (normal 0 4) 1 initial-labels) :G (dirichlet-process 10 (normal 0 4)))
(marginalfoo (normal :mu 3) (mixture [(normal 3 10) (normal (-3 2))] [0.1 0.9]))
(posterior [1 2 3] (normal :mu 3) (mixture [(normal 0 3) (normal 2 1)] [0.5 0.5]))

(defn predictive [c observations concentration base-measure]
  (let [crp (cr-process c concentration 0)
        newc (sample (sample crp ))]
    (if (contains? (:colours crp) newc)
      (let [idx (positions #(== newc %) c)
            yc (for [j idx] (nth observations j))
            posteriorc (posterior yc base-measure)]
        (sample posteriorc))
      (sample base-measure)
      )))

(def proposal (mixture [(normal -4 1) (normal 4 1)] [(/ 1 4) (/ 3 4)]))
(ic/view (plot/function-plot (pdf proposal) -5 5))
(ic/view (plot/histogram y))
(def initial-labels (into [] (range 0 (count y))))
(def y (concat (sample (normal -4 1) 50) (sample (normal 4 1) 50)))
(def initial-labels  (into [] (concat (into [] (take 50 (repeat 1))) (into [] (take 50 (repeat 2))))))
(def y (sample proposal 100))
(require '[incanter.charts :as plot])
(require '[incanter.core :as ic])
(ic/view (plot/function-plot (pdf proposal) -5 5))
(ic/view (plot/histogram y :density true :nbins 20))

(defn rn-alg3 [observations likelihood base-measure concentration niter]
  (let [n (count observations)
        initial-label (into [] (range 0 n))
        labels (take niter (iterate (partial update-labels observations likelihood base-measure concentration) initial-labels))]))
(update-labels y (normal 0 1) (normal 0 4) 1 initial-labels)
(def initial-labels  (into [] (take 100 (repeat 1))))
(update-labels)
(def c-draws (take 100 (iterate (partial update-c (into [] y) (normal 0 1) (normal 0 4)) initial-labels)))
(last c-draws)
(map (comp count distinct) c-draws)

(ic/view (plot/histogram no :nbins 20))
(ic/view (plot/histogram (map predictive c-draws (repeat y) (repeat concentration) (repeat (normal 0 4)))))




(def no (map read-string (clojure.string/split-lines (slurp "no"))))
(def yes (map read-string (clojure.string/split-lines (slurp "yes"))))
(def likelihood (poisson :rate))
(def base-measure (gamma 2 0.1))
(def concentration 1)
(def initial-labels (into [] (range 0 (count yes))))
(def initial-labels (into [] (take (count yes) (repeat 1))))
(def observations yes)
(update-labels observations likelihood base-measure concentration initial-labels)
(def label-draws (take 1000 (iterate (partial update-labels observations likelihood base-measure concentration) initial-labels)))
(def param-draws (map (partial update-params observations likelihood base-measure concentration ) label-draws))
(def pred-fn-draws (map (fn [params] (pdf (marginal likelihood (posterior-predictive params :G (dirichlet-process concentration base-measure))))) param-draws))
(def seed (plot/time-series-plot (range 0 70) (map (first pred-fn-draws) (range 0 50))))
(def allplots (reduce (fn [chart newfun] (plot/add-lines chart (range 0 70) (map newfun (range 0 70)))) seed (rest pred-fn-draws)))
(ic/view allplots)
(mean (marginal likelihood base-measure))




(sample (posterior-predictive (sample (poisson 10) 10) :G (dirichlet-process 1 (normal 0 1))))
(posterior-predictive (sample (poisson 10) 100) :G (dirichlet-process 3 (normal 0 1)))
(marginal :G (dirichlet-process 10 (normal 0 1)))
(posterior [1.1 2.2 2.2 3.3 3.3 3.3] :G (dirichlet-process 1 (normal 0 1)))
(def freq (frequencies [1.1 2.2 2.2 2.2 3.3 3.3]))
(discrete-real (keys freq) (vals freq))
(ic/view (plot/histogram (map (fn [x] (sample (poisson x))) (sample base-measure 10000)) :nbins 120 :density true))

(ic/view (plot/scatter-plot (range 0 70) (map (pdf (marginal likelihood base-measure)) (range 0 70))))
(ic/view (plot/add-points (plot/histogram (map (fn [x] (sample (poisson x))) (sample base-measure 20000)) :nbins 140 :density true) (range 0 70) (map (pdf (marginal likelihood base-measure)) (range 0 70))))
