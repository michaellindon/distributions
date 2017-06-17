(in-ns 'distributions.core)

(defn kolmogorov-smirnov [x y]
  (.kolmogorovSmirnovTest (new org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest) (double-array x) (double-array y)))

