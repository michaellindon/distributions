(in-ns 'distributions.core)

(defn normal-laplace
  "N(x|mu, s2 gamma)Laplace(x| tau/s)"
  [mu gamma tau s2]
  (let [s (sqrt s2)
        mu+ (- mu (* gamma s tau))
        mu- (+ mu (* gamma s tau))
        v (* (sqrt gamma) s)
        v2 (square v)
        w- (/ (probit (negate (/ mu- v))) (pdf (normal mu- v2) 0))
        w+ (/ (probit (/ mu+ v)) (pdf (normal mu+ v2) 0))
        w  (/ w- (+ w- w+))]
    (mixture [(truncated (normal mu- v2) Double/NEGATIVE_INFINITY 0) (truncated (normal mu+ v2) 0 Double/POSITIVE_INFINITY)] [w (- 1 w)])))
