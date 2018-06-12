;;;; package.lisp

(defpackage #:rf-utilities
  (:nicknames #:rf-utils)
  (:use #:cl)
  (:export #:sign
           #:vector-linear-combination
	   #:matrix-row
	   #:matrix-column
           #:linearize-matrix
	   #:argument
           #:angles-unitary))

(defpackage #:dopri
  (:nicknames #:rf-dopri)
  (:use #:cl #:rf-utilities))

(defpackage #:rf
  (:use #:cl #:rf-utilities #:dopri))
