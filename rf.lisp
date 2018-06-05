;;;; rf.lisp

(in-package #:rf)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;;;;; ;;
;; ;; Method Application ;; ;;
;; ;;;;;;;;;;;;;;;;;;;;;;;; ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;
;; General ;;
;;;;;;;;;;;;;

(defparameter *n-equations* 2 "Number of ODEs in the system")
(defparameter *mc* 1.0d0)
(defparameter *charge-of-particle* -1.0d0
  "Charge of particle. -1 for electron.")

(defun monodromy (n-repetition e-amplitude sigma n-cycles cep time-delay k-perpendicular ikmax kparmin kparmax
                  &key (basic-shape 11) (epsabs0 1.0e-13) (stream *standard-output*))
  "Driver function for the monodromy method."
  (loop
     with xin = (* -0.5d0 n-repetition time-delay)
     with xend = (- xin)
     with epsilonabs = epsabs0
     with cep = (* cep pi)
     for ik from 0 to ikmax
     for k-parallel = (+ kparmin (* (- kparmax kparmin) (float ik) (/ ikmax)))
     do
       (multiple-value-bind (monodromy theta1 theta2 gamma beta eigenvectors iflag)
           (monodromy-matrix xin xend epsilonabs e-amplitude sigma n-cycles k-perpendicular cep basic-shape)
         (let* ((prob-particle (expt (abs (aref monodromy 0 0)) 2))
                (prob-antiparticle (expt (abs (aref monodromy 1 0)) 2))
                (unit (abs (- 1.0d0 prob-particle prob-antiparticle)))
                (sin2theta (max (+ (expt (abs (aref monodromy 0 1)) 2)
                                   (* 0.25 (expt (abs (- (aref monodromy 0 0)
                                                         (aref monodromy 1 1))) 2)))
                                1.0d-5))
                (sin2gamma (/ (expt (abs (aref monodromy 0 1)) 2) sin2theta)))
           (format stream "~20,14e ~20,14e ~20,14e ~20,14e ~20,14e ~20,14e ~20,14e ~20,14e ~20,14e ~20,14e~
                           ~d ~{~20,14e~^ ~} ~{~20,14e~^ ~}" k-parallel prob-particle prob-antiparticle sin2gamma sin2theta
                           theta1 theta2 gamma beta unit
                           iflag (linearize-matrix monodromy :unpack-complex t)
                           (linearize-matrix eigenvectors :unpack-complex t))))))

(defun monodromy-matrix (xin xend epsilonabs e-amplitude sigma n-cycles k-perpendicular cep basic-shape)
  (let* ((monodromy (make-array `(,*n-equations* ,*n-equations*) :initial-element 0))
         (sout (make-array *n-equations* :initial-element 0))
         ;; (yend (make-array (* 2 *n-equations*) :initial-element 0))
         (eigenvectors (make-array `(,*n-equations* ,*n-equations*) :initial-element 0))
         (ysystem (make-array (* 2 *n-equations*) :initial-element 0))
         (ypsystem (make-array (* 2 *n-equations*) :initial-element 0))
         (jumpmax (/ (- xend xin) 1000.0))
         (jumpguess (/ jumpmax 100.0))
         (nscale (/ (- xend xin) jumpguess))
         (scaling-factor (scaling-factor xin xend nscale basic-shape sigma n-cycles cep)))
    (if (< (abs scaling-factor) (* 100.0 double-float-epsilon))
        (error "Incorrect field parameters. Scaling factor too small (value: ~a)" scaling-factor)
        (loop
           for column from 0 below *n-equations*
           for yend = (make-array (* 2 *n-equations*) :initial-element 0)
           do
             (setf (aref yend (* 2 column)) 1.0)))
    (values monodromy theta1 theta2 gamma beta eigenvectors iflag)))


(defun scaling-factor (tscalemin tscalemax nscale basic-shape sigma n-cycles cep)
  (loop
     with scaling-factor = 0.0
     for is from 0 to nscale
     for tt = (+ tscalemin (* 1.0 is (/ (- tscalemax tscalemin) nscale)))
     for value = (abs (electric0 basic-shape tt sigma n-cycles cep))
     do (when (> value scaling-factor)  ; TODO: Replace with MAX
          (setf scaling-factor value))
     finally (return scaling-factor)))

;; Generic function approach to calculating the electric field vector

(defgeneric electric0 (shape tt sigma n-cycles cep))

(defmethod electric0 ((shape (eql 11)) tt sigma n-cycles cep)
  (if (< (abs tt) sigma)
      (+ (* 0.5 (sin (+ (/ (* 2.0 n-cycles pi tt) sigma) cep)))
         (* 0.25 (sin (+ (/ (* (+ (* 2.0 n-cycles) 1.0) pi tt) sigma) cep)))
         (* 0.25 (sin (+ (/ (* (+ (* 2.0 n-cycles) -1.0) pi tt) sigma) cep))))
      0.0))

;;;;;;;;;;;;;;;;;;;;;;;;;
;; Starting parameters ;;
;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *i20rf1* '(1              ;n-repetition
                         0.1            ;e-amplitude
                         10.0           ;sigma
                         2              ;n-cycles
                         0.0            ;cep
                         200.0          ;time-delay
                         0.0            ;k-perpendicular
                         10000          ;ikmax
                         -2.0           ;kparmin
                         2.0            ;kparmax
                         :epsabs0 1.0e-13))

(defparameter *i20rf2* '(2              ;n-repetition
                         0.1            ;e-amplitude
                         10.0           ;sigma
                         2              ;n-cycles
                         0.0            ;cep
                         200.0          ;time-delay
                         0.0            ;k-perpendicular
                         10000          ;ikmax
                         -2.0           ;kparmin
                         2.0            ;kparmax
                         :epsabs0 1.0e-13))

(defparameter *i20rf3* '(3              ;n-repetition
                         0.1            ;e-amplitude
                         10.0           ;sigma
                         2              ;n-cycles
                         0.0            ;cep
                         200.0          ;time-delay
                         0.0            ;k-perpendicular
                         10000          ;ikmax
                         -2.0           ;kparmin
                         2.0            ;kparmax
                         :epsabs0 1.0e-13))
