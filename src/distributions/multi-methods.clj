(in-ns 'distributions.core)

(defmulti posterior (fn [data likelihood prior] [(class likelihood) (class prior)]))

(defmulti marginal (fn [likelihood prior] [(class likelihood) (class prior)]))




