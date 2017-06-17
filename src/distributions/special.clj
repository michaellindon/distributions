(in-ns 'distributions.core)
(import '(org.apache.commons.math3.special Beta Gamma Erf))

(defn gamma-fn [x] (Gamma/gamma x))
(defn digamma-fn [x] (Gamma/digamma x))
(defn trigamma-fn [x] (Gamma/trigamma x))
(defn log-gamma-fn [x] (Gamma/logGamma x))
(defn log-beta-fn [a b] (Beta/logBeta a b))
(defn beta-fn [a b] (Math/exp (Beta/logBeta a b)))
(defn erf
  ([x] (Erf/erf x))
  ([a b] (Erf/erf a b)))
(defn erfc [x] (Erf/erfc x))
(defn inverse-erf [x] (Erf/erfInv x))
(defn inverse-erfc [x] (Erf/erfcInv x))
