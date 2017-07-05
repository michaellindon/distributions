(in-ns 'distributions.core)

(defn gamma-fn [x] (org.apache.commons.math3.special.Gamma/gamma x))
(defn digamma-fn [x] (org.apache.commons.math3.special.Gamma/digamma x))
(defn trigamma-fn [x] (org.apache.commons.math3.special.Gamma/trigamma x))
(defn log-gamma-fn [x] (org.apache.commons.math3.special.Gamma/logGamma x))
(defn log-beta-fn [a b] (org.apache.commons.math3.special.Beta/logBeta a b))
(defn beta-fn [a b] (Math/exp (org.apache.commons.math3.special.Beta/logBeta a b)))
(defn erf
  ([x] (org.apache.commons.math3.special.Erf/erf x))
  ([a b] (org.apache.commons.math3.special.Erf/erf a b)))
(defn erfc [x] (org.apache.commons.math3.special.Erf/erfc x))
(defn inverse-erf [x] (org.apache.commons.math3.special.Erf/erfInv x))
(defn inverse-erfc [x] (org.apache.commons.math3.special.Erf/erfcInv x))
