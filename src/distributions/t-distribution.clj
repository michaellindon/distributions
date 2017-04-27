(in-ns 'distributions.core)

(defn t-distribution
  ([df]
   (new TDistribution df))
  ([df location scale]
   (location-scale (t-distribution df) location scale)))

