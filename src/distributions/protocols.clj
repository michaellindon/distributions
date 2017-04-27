(in-ns 'distributions.core)

(defprotocol probability-function
  (P [d] [d x]))

(defprotocol density-function
  (pdf [d] [d x])
  (log-pdf [d] [d x]))

(defprotocol mass-function
  (pmf [d] [d x])
  (log-pmf [d] [d x]))

(defprotocol distribution-function
  (cdf [d] [d x]))

(defprotocol support
  (support-lower [d])
  (support-upper [d]))

(defprotocol inverse-distribution-function
  (icdf [d] [d x]))

(defprotocol first-moment
  (mean [d]))

(defprotocol second-central-moment
  (variance [d]))

(defprotocol standard-deviation
  (std [d]))

(defprotocol proximal
  (prox [d] [d h x]))

(defprotocol random
  (sample [d] [d n]))
