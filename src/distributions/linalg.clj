(in-ns 'distributions.core)

(defn log-det [M]
  (let [L (:L (la/cholesky M))]
    (* 2 (reduce + (map log (diagonal L))))))

