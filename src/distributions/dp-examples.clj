(ns distributions.dp-examples
  (:require [distributions.core :refer :all]))


;Normal Mixture Model Example
(require '[incanter.charts :as plot])
(require '[incanter.core :as ic])
(def proposal (mixture [(normal -4 1) (normal 4 1)] [(/ 1 4) (/ 3 4)]))
(def y (sample proposal 100))
(def initial-labels (into [] (range 0 (count y))))
(def likelihood (normal :mu 1))
(def base-measure (normal 0 4))
(def observations y)
(def concentration 0.1)
(def label-draws (take 100 (iterate (partial dp-update-labels observations likelihood base-measure concentration) initial-labels)))
(def param-draws (map (partial dp-update-params observations likelihood base-measure concentration ) label-draws))
(def pred-fn-draws (map (fn [params] (pdf (marginal (normal :mu 1) (posterior-predictive params :G (dirichlet-process concentration base-measure))))) param-draws))
(def seed (plot/function-plot (first pred-fn-draws) -5 5))
(def allplots (reduce (fn [chart newfun] (plot/add-function chart newfun -5 5)) seed (rest pred-fn-draws)))
(ic/view allplots)
(rn-alg3 observations likelihood base-measure concentration 10)








;Traffic Examples
(def no (map read-string (clojure.string/split-lines (slurp "no"))))
(def yes (map read-string (clojure.string/split-lines (slurp "yes"))))
(def likelihood (poisson :rate))
(def base-measure (gamma 2 0.1))
(def concentration 1)
(def initial-labels (into [] (range 0 (count yes))))
(def observations yes)
(def label-draws (take 100 (iterate (partial dp-update-labels observations likelihood base-measure concentration) initial-labels)))
(def param-draws (map (partial dp-update-params observations likelihood base-measure concentration ) label-draws))
(def pred-fn-draws (map (fn [params] (pdf (marginal likelihood (posterior-predictive params :G (dirichlet-process concentration base-measure))))) param-draws))
(def seed (plot/time-series-plot (range 0 70) (map (first pred-fn-draws) (range 0 50))))
(def allplots (reduce (fn [chart newfun] (plot/add-lines chart (range 0 70) (map newfun (range 0 70)))) seed (rest pred-fn-draws)))
(ic/view allplots)




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
(sample (posterior-predictive (sample (poisson 10) 10) :G (dirichlet-process 1 (normal 0 1))))
(posterior-predictive (sample (poisson 10) 100) :G (dirichlet-process 3 (normal 0 1)))
(marginal :G (dirichlet-process 10 (normal 0 1)))
(posterior [1.1 2.2 2.2 3.3 3.3 3.3] :G (dirichlet-process 1 (normal 0 1)))
(def freq (frequencies [1.1 2.2 2.2 2.2 3.3 3.3]))
(discrete-real (keys freq) (vals freq))
(ic/view (plot/histogram (map (fn [x] (sample (poisson x))) (sample base-measure 10000)) :nbins 120 :density true))

(ic/view (plot/scatter-plot (range 0 70) (map (pdf (marginal likelihood base-measure)) (range 0 70))))
(ic/view (plot/add-points (plot/histogram (map (fn [x] (sample (poisson x))) (sample base-measure 20000)) :nbins 140 :density true) (range 0 70) (map (pdf (marginal likelihood base-measure)) (range 0 70))))





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
