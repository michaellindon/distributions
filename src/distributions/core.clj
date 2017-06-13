(ns distributions.core
  (:require [distributions.root :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.core.matrix.linear :as la]))

(set-current-implementation :vectorz)
(defn inv [x] (/ 1 x))
(load "protocols")
(load "acm-distributions")
(def probit (cdf (normal 0 1)))
(load "inverse-gaussian")
(load "gdp")
(load "defaults")
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

