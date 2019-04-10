;;;; rf-test.lisp

(in-package #:rf)

;;;;;;;;;;;;;;;;;;
;; Test Drivers ;;
;;;;;;;;;;;;;;;;;;

(defun test-harmonic-oscillator-bld (&key (stream *standard-output*))
  ""
  (let* ((xin 0.0d0)
         (yin (make-array '(2) :initial-contents
                          (list 1.0d0 0.0d0)))
         (xend 10.0d0)
         (epsilonabs 1.0d-15)
         (jumpmax (/ (- xend xin) 100.0d0))
         (jumpguess (/ jumpmax 10.0d0))
         (n-equations 1)
         (fn-misc nil)
         (trail 
          (bld-ode::rka #'f-harmonic-oscillator xin xend yin
                        :hmax jumpmax :h0 jumpguess
                        :tol epsilonabs))
         (yend (cadar (last trail))))
    (format stream "~20,14e ~20,14e~%~
~20,14e ~20,14e~%~
~20,14e ~20,14e~%"
            (aref yend 0) (aref yend 1)
            (cos xend) (sin xend)
            (abs (- (aref yend 0) (cos xend)))
            (abs (- (aref yend 1) (sin xend))))))

(defun test-harmonic-oscillator-driver (&key (stream *standard-output*))
  ""
  (let* ((xin 0.0d0)
         (yin (make-array '(2) :initial-contents
                          (list 1.0d0 0.0d0)))
         (xend 10.0d0)
         (epsilonabs 1.0d-15)
         (jumpmax (/ (- xend xin) 100.0d0))
         (jumpguess (/ jumpmax 10.0d0))
         (n-equations 1)
         (fn-misc nil)
         (yend (dopri::dopri8
                (* 2 n-equations) #'f-harmonic-oscillator
                fn-misc xin yin xend epsilonabs jumpmax jumpguess)))
    (format stream "~20,14e ~20,14e~%~
~20,14e ~20,14e~%~
~20,14e ~20,14e~%"
            (aref yend 0) (aref yend 1)
            (cos xend) (sin xend)
            (abs (- (aref yend 0) (cos xend)))
            (abs (- (aref yend 1) (sin xend))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; System-definition functions ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun f-harmonic-oscillator (tt y misc)
  "ODE system for harmonic oscillator
   
   TT - independent variable
   Y - N-dimensional initial vector for dependent variable"
  (declare (ignorable tt misc))
  (declare (type (simple-vector *) y))
  (let* ((n (array-total-size y))
         (yend (make-array `(,n) :initial-element 0.0d0)))
    (setf (aref yend 0)
          (* -1.0d0 (aref y 1)))
    (setf (aref yend 1)
          (* 1.0d0 (aref y 0)))
    yend))
