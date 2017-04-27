(ns distributions.root)

(defn bisection
  ([f a b]
   (bisection f a b (* 100 (java.lang.Math/ulp 1.0))))
  ([f a b eps]
   (loop [l a
          u b
          fl (f l)
          fu (f u)]
     (let [m (* 0.5 (+ l u))
           fm (f m)]
       (if (< (Math/abs fm) eps)
         m
         (if (= (java.lang.Math/signum fm) (java.lang.Math/signum fl)) (recur m u fm fu) (recur l m fl fm)))))))

(defn newton-raphson-step [f f' x] (- x (/ (f x) (f' x))))

(defn newton-raphson
  ([f f' initial-x]
   (newton-raphson f f' initial-x (java.lang.Math/ulp 1.0)))
  ([f f' initial-x eps]
   (let [nr-iterator (partial newton-raphson-step f f')
         nr-stream (iterate nr-iterator initial-x)]
     (first (drop-while #(> (Math/abs (f %)) eps) nr-stream))))
  )

(defn secant-step [f [x-1 x-2]] [(- x-1 (* (f x-1) (/ (- x-1 x-2) (- (f x-1) (f x-2))))) x-1])

(defn secant
  ([f x0 x1]
   (secant f x0 x1 (java.lang.Math/ulp 1.0)))
  ([f x0 x1 eps]
   (let [g (memoize f)
         s-iterator (partial secant-step g)
         s-stream (iterate s-iterator [x1 x0])]
     (first (first (drop-while #(> (Math/abs (f (first %))) eps) s-stream))))))
