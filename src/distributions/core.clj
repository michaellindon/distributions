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

(comment 
(gg4clj/render [[:<- :g (gg4clj/data-frame g-dat)]
                (gg4clj/r+
                 [:ggplot :g [:aes :g1 :g2]]
                 [:xlim -2 2]
                 [:ylim -2 2]
                 [:geom_point {:colour "steelblue" :size 2}]
                 [:stat_density2d {:colour "#FF29D2"}]
                 [:theme_bw])]
               {:width 5 :height 5})

(gg4clj/render [[:<- :g (gg4clj/data-frame g-dat)]
                (gg4clj/r+
                 [:ggplot :g [:aes :g1]]
                 [:geom_density {:adjust 1 :colour "steelblue"}]
                 [:theme_bw]
                 [:stat_function {:function exp}]
                 )]
               {:width 5 :height 5})

(def foo (sample (normal 0 1) 100))
(gg4clj/render (gg4clj/to-r [:qplot foo]))
(gg4clj/render [[:<- :g (gg4clj/data-frame g-dat)]
              (gg4clj/r+
               [:stat_function {:fun (pdf the-d) :colour "red"}]
               )]
             {:width 5 :height 5})

(def foo (new org.rosuda.REngine.Rserve.RConnection))

(defn randJ [z]
  (let [t 0.64
        ex (exponential (+ (* 0.5 (square z)) (/ (square Math/PI) 8)))
        ig (inverse-gaussian (inv (abs z)) 1)
        p (* 2 (cosh z) (exp (negate z)) (cdf ig t))
        q (* 0.5 Math/PI (- 1 (cdf ex t)))
        ig-prob (/ p (+ p q))
        cg (fn [x] (*
                    (* 0.5 Math/PI)
                    (cosh z)
                    (if (< x t)
                      (* (pow (/ 2 (* x Math/PI)) 1.5) (exp (negate (* 0.5 (+ (/ 1 x) (* x (square z))))) ))
                      (exp (+ (* 0.5 (square z)) (/ (square Math/PI) 8))))))
        ]
    (loop []
      (let [X (sample (mixture [(truncated ig 0 t) (truncated ex t Double/POSITIVE_INFINITY)] [ig-prob (- 1 ig-prob)]))
            U (sample (uniform 0 (cg X)))
            stream (iterate (partial jacobi-iterator X z) [(jacobi-a X 0 z) 0])
            S (first (drop-while (fn [[sn n]] (if (odd? n) (> U sn) (< U sn))) stream))]
        ;(if (odd? (second S)) X (recur))
        X
        )
      )
    ))

(truncated (inverse-gaussian 1 1) 0 0.64)

(defn jacobi-a
  ([x n z] (* (cosh z) (exp (negate (* z z 0.5 x))) (jacobi-a x n)))
  ([x n]
   (if (< x 0.64)
     (* Math/PI (+ n 0.5) (pow (/ 2 (* x Math/PI)) 1.5) (exp (negate (/ (* 2 (square (+ n 0.5))) x))))
     (* Math/PI (+ n 0.5) (exp (negate (* 0.5 x (square (* Math/PI (+ n 0.5)))))) )))
  )

(defn not-converged? [[[sn _] [sp _]]] (> (abs (- sn sp)) (java.lang.Math/ulp 1.0)))

(def t 0.64)
(defn jacobi-iterator [x z [s n]] [(if (even? n) (+ s (jacobi-a x (inc n) z)) (- s (jacobi-a x (inc n) z)) ) (inc n) z])
(defn jacobi-density [z x]
  (let [stream (iterate (partial jacobi-iterator x z) [(jacobi-a x 0 z) 0])
        zip-stream (map vector (drop 1 stream) stream)]
    (first (first (first (drop-while not-converged? zip-stream))))))

(randJ 1)
(def z 0.01)
(defn jacobi-density [x] (first (last (take 100 (iterate jacobi-iterator [(jacobi-a x 0 z) 0 z])))))
(take 3 (iterate jacobi-iterator [(jacobi-a t x 0) 0]))
(jacobi-density 1)
(qplot (range 0.01 4 0.01) (map (fn [x] (jacobi-density x 1.5)) (range 0.01 4 0.01)))
(close-plot 2)

(jacobi-a 0.64 0.3 4)
(def t 0.64)
(mixture [])


(.eval foo "library(ggplot2)")
(.eval foo "print(qplot(x=rnorm(100),y=rnorm(100)))")
(.eval foo "x=print(ggplot(diamonds, aes(depth, fill = cut, colour = cut))+geom_density(alpha = 0.1))")
(.eval foo "x=qplot(x=rnorm(100),y=rnorm(100))")
(.eval foo "x+geom_density(alpha = 0.1)")
(.eval foo "print(x + xlim(-200,200)+stat_function(fun=dnorm))")
(.eval foo "dev.off()")
(.eval foo "plot(x=rnorm(100),y=rnorm(100))")
(.voidEval foo "bob <- 20")
(.eval foo "ban")
(.eval foo "{catty=ggplot(diamonds, aes(carat)) +
  geom_density(); catty}")

(take 2000 (repeatedly #(randJ 1.5)))

(.assign foo "y" "bar")
(.eval foo "x <- 1001")
(.asString (.eval foo "y"))
(.assign foo "bar" "e")
(.get foo "bar")
(def bar (double-array [1 2 3]))
(.eval foo "plot(y)")
(.eval foo "y <- 100*x")
(vec (.asDoubles (.eval foo "y")))
(.assign foo "x" (double-array (take 2000 (repeatedly #(randJ 1.5))) ))
(.parseAndEval foo "print(qplot(x,geom=\"histogram\"))")
(.toString (.asNativeJavaObject (first (.eval foo "print(x)"))))
(defn qplot [x y]
  (.eval foo "dev.new()")
  (.assign foo "x" (double-array x))
  (.assign foo "y" (double-array y))
  (.eval foo "print(qplot(x,y))")
  )

(defn switch-plot [x]
  (.eval foo (str "dev.set(" x ")")))
(defn close-plot [x]
  (.eval foo (str "dev.off(" x ")")))

(cdf (inverse-gaussian 3 2) -0.1)
(qplot (range 0 100) (map sin (range 0 100)))
(switch-plot 3)


(close-plot 3)
(.eval foo "dev.off()")

)
