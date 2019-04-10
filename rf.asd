;;;; rf.asd

(asdf:defsystem #:rf
  :description "Implementation of the Dormand-Prince algorithm for
  solving nonstiff systems of differential equations."
  :author "Wojciech Gac <wojciech.s.gac@gmail.com>"
  :license  "GPLv3"
  :version "0.0.1"
  :serial t
  :depends-on (;; #:matplotlib-cl
               #:computable-reals
               #:bld-ode)
  :components ((:file "package")
               (:file "utilities")
               (:file "dopri")
               (:file "rf")
               (:file "rf-test")))
