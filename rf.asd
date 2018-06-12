;;;; rf.asd

(asdf:defsystem #:rf
  :description "Describe rf here"
  :author "Wojciech Gac <wojciech.s.gac@gmail.com>"
  :license  "GPLv3"
  :version "0.0.1"
  :serial t
  :depends-on (;; #:matplotlib-cl
               #:computable-reals)
  :components ((:file "package")
               (:file "utilities")
               (:file "dopri")
               (:file "rf")))
