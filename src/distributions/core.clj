(ns distributions.core
  (:require [distributions.root :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.linear :as la]
            [clojure.set :as set]))

(set-current-implementation :vectorz)

(defn normalize [coll]
  (if (empty? coll)
    nil
    (let [Z (reduce + coll)]
      (map #(/ % Z) coll))))

(defn normalize-log [coll]
  (let [a (reduce max coll)
        shifted (map - coll (repeat a))
        expshifted (map exp shifted)
        Z (reduce + expshifted)
        ]
    (map #(/ % Z) expshifted)))

(defn log-sum-exp [coll]
  (let [a (reduce max coll)
        expshifted (map (fn [x] (exp (- x a))) coll)]
    (+ a (log (reduce + expshifted)))))

(defn inv [x] (/ 1 x))
(defn positions
  [pred coll]
  (keep-indexed (fn [idx x]
                  (when (pred x)
                    idx))
                coll))

(defn remove-at [coll idx] (for [i (remove #(== % idx) (range 0 (count coll)))] (nth coll i)))
(load "multi-methods")
(load "linalg")
(load "protocols")
(load "inverse-gamma")
(load "normal")
(def probit (cdf (normal 0 1)))
(load "gamma")
(load "negative-binomial")
(load "poisson")
(load "acm-distributions")
(load "inverse-gaussian")
(load "gdp")
(load "defaults")
(load "discrete-integer")
(load "mixture")
(load "mvnormal")
(load "truncated")
(load "location-scale")
(load "t-distribution")
(load "normal-laplace")
(load "polya-gamma")
(load "integration")
(load "special")
(load "dirichlet")
(load "discrete-real")
(load "mvt")
(load "inference")
(load "enumerated")
(load "chinese-restaurant-process")

(defn posterior-predictive [data likelihood prior]
  (marginal likelihood (posterior data likelihood prior)))

(load "dirichlet-process")
(load "sampling")
(load "normal-inverse-gamma")
