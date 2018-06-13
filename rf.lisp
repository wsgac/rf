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

(defun monodromy-driver (n-repetition e-amplitude sigma n-cycles cep time-delay k-perpendicular ikmax kparmin kparmax
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
           (monodromy-matrix xin xend epsilonabs e-amplitude sigma n-cycles n-repetition
                             k-perpendicular k-parallel time-delay cep basic-shape)
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

(defun monodromy-matrix (xin xend epsilonabs e-amplitude sigma n-cycles n-repetition
                         k-perpendicular k-parallel time-delay cep basic-shape)
  (let* ((monodromy (make-array `(,*n-equations* ,*n-equations*) :initial-element 0))
         ;; (sout (make-array *n-equations* :initial-element 0))
         ;; (yend (make-array (* 2 *n-equations*) :initial-element 0))
         ;; (eigenvectors (make-array `(,*n-equations* ,*n-equations*) :initial-element 0))
         ;; (ysystem (make-array (* 2 *n-equations*) :initial-element 0))
         ;; (ypsystem (make-array (* 2 *n-equations*) :initial-element 0))
         (jumpmax (/ (- xend xin) 1000.0))
         (jumpguess (/ jumpmax 100.0))
         (nscale (/ (- xend xin) jumpguess))
         (scaling-factor (scaling-factor xin xend nscale basic-shape sigma n-cycles cep))
         (fn-misc (list
                   ;; Parameter packing for F-SYSTEM-REAL
                   :shape basic-shape
                   :k-perpendicular k-perpendicular
                   :k-parallel k-parallel
                   :e-amplitude e-amplitude
                   :scaling-factor scaling-factor
                   :time-delay time-delay
                   :n-repetition n-repetition
                   :n-cycles n-cycles
                   :sigma sigma
                   :cep cep)))
    (if (< (abs scaling-factor) (* 100.0 double-float-epsilon))
        (error "Incorrect field parameters. Scaling factor too small (value: ~a)" scaling-factor)
        (loop
           for column from 0 below *n-equations*
           for yend = (make-array (* 2 *n-equations*) :initial-element 0)
           do
             (setf (aref yend (* 2 column)) 1.0)
             (let ((yend (dopri::dopri8
                             (* 2 *n-equations*) #'f-system-real fn-misc xin yend xend epsilonabs jumpmax jumpmax)))
               (loop
                  for ic from 0 below *n-equations*
                  do
                    ;; (setf (aref sout ic)
                    ;;       (complex (aref yend (* 2 ic)) (aref yend (1+ (* 2 ic)))))
                    (setf (aref monodromy ic column)
                          (complex (aref yend (* 2 ic)) (aref yend (1+ (* 2 ic)))))))
           finally
             (let ((a (aref monodromy 0 0))
                   (b (aref monodromy 0 1))
                   (c (aref monodromy 1 0))
                   (d (aref monodromy 1 1))
                   (eps 1.0d-7))
               (multiple-value-bind (theta1 theta2 gamma beta eigenvectors)
                   (angles-unitary a b c d eps)
                 (values monodromy theta1 theta2 gamma beta eigenvectors ;; iflag
                         )))))))


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

(defun electric (shape tt e-amplitude sigma time-delay n-cycles n-repetition scaling-factor cep)
  (loop
     for i from 1 to n-repetition
     sum (electric0 shape
                    (+ tt (* (- i (* 0.5d0 (+ n-repetition 1.0d0))) time-delay))
                    sigma
                    n-cycles
                    cep) into electric
     finally (return (* e-amplitude electric (/ scaling-factor)))))

;; Vector potential function - A(t)
(defun a-field (shape tt e-amplitude scaling-factor time-delay n-repetition n-cycles sigma cep)
  (loop
     for i from 1 to n-repetition
     sum (a-field-0 shape
                    (+ tt (* (- i (* 0.5d0 (+ n-repetition 1.0d0))) time-delay))
                    sigma
                    n-cycles
                    cep)
     into afield
     finally (return (* e-amplitude afield (/ scaling-factor)))))

(defgeneric a-field-0 (shape tt sigma n-cycles cep))

(defmethod a-field-0 ((shape (eql 11)) tt sigma n-cycles cep)
  (if (>= (abs tt) sigma)
      0.0d0
      (+ (* 0.25d0 (/ sigma (* pi n-cycles))
            (- (cos (+ (* 2.0d0 n-cycles pi (/ tt sigma)) cep))
               (cos cep)))
         (* 0.125d0 (/ sigma (* pi (+ n-cycles 1.0d0)))
            (+ (cos (+ (* (+ (* 2.0d0 n-cycles) 1.0d0) pi (/ tt sigma)) cep))
               (cos cep)))
         (* 0.125d0 (/ sigma (* pi (- n-cycles 1.0d0)))
            (+ (cos (+ (* (- (* 2.0d0 n-cycles) 1.0d0) pi (/ tt sigma)) cep))
               (cos cep))))))

;; Helper functions called when integrating

(defun f-system-real (n tt y misc)
  "Helper function for SYSTEM-REAL. Its main purpose is to unpack
   various assorted parameters from MISC and pass them to
   SYSTEM-REAL."
  (declare (ignorable n))
  (let ((shape (getf misc :shape))
        (k-perpendicular (getf misc :k-perpendicular))
        (k-parallel (getf misc :k-parallel))
        (e-amplitude (getf misc :e-amplitude))
        (scaling-factor (getf misc :scaling-factor))
        (time-delay (getf misc :time-delay))
        (n-repetition (getf misc :n-repetition))
        (n-cycles (getf misc :n-cycles))
        (sigma (getf misc :sigma))
        (cep (getf misc :cep)))
    (system-real shape tt y k-perpendicular k-parallel e-amplitude scaling-factor
                 time-delay n-repetition n-cycles sigma cep)))

(defun system-real (shape tt y-in k-perpendicular k-parallel e-amplitude scaling-factor
                    time-delay n-repetition n-cycles sigma cep)
  "Actual function defining the system of ODEs."
  (let ((y-out (make-array (* 2 *n-equations*) :initial-element 0))
        (f-epsilon (sqrt (func-epsilon-2 shape tt k-perpendicular k-parallel e-amplitude
                                         scaling-factor time-delay n-repetition n-cycles sigma cep)))
        (f-omega (func-omega shape tt k-perpendicular k-parallel e-amplitude
                             scaling-factor time-delay n-repetition n-cycles sigma cep)))
    (setf (aref y-out 0) (+ (* f-epsilon (aref y-in 1)) (* f-omega (aref y-in 2))))
    (setf (aref y-out 1) (+ (- (* f-epsilon (aref y-in 0))) (* f-omega (aref y-in 3))))
    (setf (aref y-out 2) (+ (- (* f-epsilon (aref y-in 3))) (- (* f-omega (aref y-in 0)))))
    (setf (aref y-out 3) (+ (* f-epsilon (aref y-in 2)) (- (* f-omega (aref y-in 1)))))
    y-out))

(defun func-epsilon-2 (shape tt k-perpendicular k-parallel e-amplitude scaling-factor time-delay
                       n-repetition n-cycles sigma cep)
  ""
  (+ (expt *mc* 2) (expt k-perpendicular 2)
     (expt (- k-parallel (* *charge-of-particle*
                            (a-field shape tt e-amplitude scaling-factor time-delay
                                     n-repetition n-cycles sigma cep))) 2)))

(defun func-omega (shape tt k-perpendicular k-parallel e-amplitude scaling-factor time-delay
                   n-repetition n-cycles sigma cep)
  ""
  (* -0.5d0 *charge-of-particle* (electric shape tt e-amplitude sigma time-delay n-cycles n-repetition scaling-factor cep)
     (sqrt (+ (expt *mc* 2) (expt k-perpendicular 2)))
     (/ (func-epsilon-2 shape tt k-perpendicular k-parallel e-amplitude
                        scaling-factor time-delay n-repetition n-cycles sigma cep))))

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

;; Test

(defparameter *test* '(1              ;n-repetition
                       0.1            ;e-amplitude
                       10.0           ;sigma
                       2              ;n-cycles
                       0.0            ;cep
                       200.0          ;time-delay
                       0.0            ;k-perpendicular
                       1000          ;ikmax
                       -2.0           ;kparmin
                       2.0            ;kparmax
                       :epsabs0 1.0e-6))

(defparameter *dopri8-test*
  (list
   4
   #'f-system-real
   (list :SHAPE 11 :K-PERPENDICULAR 0.0 :K-PARALLEL -2.0 :E-AMPLITUDE 0.1 :SCALING-FACTOR 0.9630919418921946d0 :TIME-DELAY 200.0 :N-REPETITION 1 :N-CYCLES 2 :SIGMA 10.0 :CEP 0.0d0)
   -100.0d0
   #(0 0 1.0 0)
   100.0d0
   1.e-6
   0.2d0
   0.2d0))

;;;;;;;;;;;;;;;;;;;;;;;
;; Wrapper functions ;;
;;;;;;;;;;;;;;;;;;;;;;;

(defun run-driver (driver-fn input-parameters output-file)
  (with-open-file (out output-file
		       :direction :output
		       :if-exists :supersede)
    (let ((*standard-output* out))
      (apply driver-fn input-parameters))))
