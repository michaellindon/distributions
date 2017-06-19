(in-ns 'distributions.core)


(defn log-det [M] (reduce + (map log (:S (la/svd M)))))

