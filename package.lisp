;;;; package.lisp

(defpackage #:rf-utilities
  (:nicknames #:rf-utils)
  (:use #:cl)
  (:export #:sign
           #:vector-linear-combination))

(defpackage #:rf-odedr-wanner-hairer
  (:nicknames #:rf-owh #:rf-odedr)
  (:use #:cl #:rf-utilities))

(defpackage #:rf
  (:use #:cl #:rf-utilities #:rf-odedr))
