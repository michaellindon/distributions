(ns distributions.integration)

(defn expectation
  ([h d n]
   (mean (map h (sample d n))))
  ([h d]
   (expectation h d 10000)))

(defn kullback-leibler
  ([f g n]
   (let [h (fn [x] (- (log-pdf f x) (log-pdf g x)))]
     (expectation h f n)))
  ([f g]
   (kullback-leibler f g 10000)))

(defn expectation-qi
  ([h d n]
   (let [delta (/ 1.0 n)
         grid (map (icdf d) (range delta 1 delta))]
     (mean (map h grid))))
  ([h d]
   (expectation-qi h d 10000)))

(defn kullback-leibler-qi
  ([f g n]
   (let [h (fn [x] (- (log-pdf f x) (log-pdf g x)))]
     (expectation-qi h f n)))
  ([f g]
   (kullback-leibler-qi 10000)))

(defn quantile-integrate
  ([f d n]
   (let [delta (/ 1.0 n)
         grid (map (icdf d) (range delta 1 delta))
         f_i (map f grid)]
     (/ (reduce + 0 f_i) n))))
