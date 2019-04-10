;;;; utilities.lisp

(in-package #:rf-utilities)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Fortran compatibility ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun sign (a b)
  "Implements a Fortran-like SIGN function. Returns the absolute value
   of A with the sign taken from B."
  (* (abs a) (if (minusp b) -1 1)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Matrix and vector functions ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun vector-linear-combination (coefficients vectors &key (type 'vector))
  "Calculate a linear combination of VECTORS and COEFFICIENTS in an
  element-wise fashion. Adjust TYPE to reflect whether VECTORS are of
  type 'VECTOR or 'LIST."
  (assert (member type '(vector list)))
  (assert (listp coefficients))
  (assert (listp vectors))
  (assert (= (length coefficients) (length vectors)))
  (assert (every #'numberp coefficients))
  (assert (every #'(lambda (x) (typep x 'sequence)) vectors))
  (assert (apply #'= (mapcar #'length vectors)))
  (apply #'map
         type
         (lambda (&rest vector-components)
           (loop
              for coefficient in coefficients
              for vc in vector-components
              sum (* coefficient vc)))
         vectors))

(defun matrix-row (matrix row)
  "Extract ROW from a two-dimensional MATRIX as a one-dimensional vector."
  (assert (= (array-rank matrix) 2) (matrix) "Not a two-dimensional matrix: ~a" matrix)
  (loop
     for col from 0 below (array-dimension matrix 1)
     collect (aref matrix row col) into res
     finally (return (coerce res 'vector))))

(defun matrix-column (matrix col)
  "Extract COLUMN from a two-dimensional MATRIX as a one-dimensional vector."
  (assert (= (array-rank matrix) 2) (matrix) "Not a two-dimensional matrix: ~a" matrix)
  (loop
     for row from 0 below (array-dimension matrix 0)
     collect (aref matrix row col) into res
     finally (return (coerce res 'vector))))

(defun linearize-matrix (matrix &key unpack-complex)
  "Convert a 2-dimensional array into a flat list of elements in
   row-major order. If UNPACK-COMPLEX is T, represent complex numbers
   as pairs of real and imaginary parts."
  (destructuring-bind (rows columns)
      (array-dimensions matrix)
    (loop
       for row from 0 below rows
       append (loop
                 for column from 0 below columns
                 append (let ((el (aref matrix row column)))
                          (if unpack-complex
                              (list (realpart el) (imagpart el))
                              (list el)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Complex-number-specific functions ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Crate an alias so that the argument of a complex number can be
;; obtained with ARGUMENT
(setf (fdefinition 'argument) #'phase)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Unitary operator functions ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun angles-unitary (a b c d eps)
  "Given the elements of a unitary matrix - A, B, C and D - determine
   the angles THETA1 and THETA2 as well as EIGENVECTORS and
   coefficients BETA and GAMMA. Return them as multiple values."
  (let ((flag 0)
	(eigenvectors (make-array `(,rf::*n-equations* ,rf::*n-equations*)
				  :initial-element 0)))

    (when (> (abs (- (+ (expt (abs a) 2) (expt (abs c) 2)) 1.0)) eps)
      (incf flag 10))
    (when (> (abs (- (+ (expt (abs b) 2) (expt (abs d) 2)) 1.0)) eps)
      (incf flag 100))
    (when (> (abs (+ (* b (conjugate a)) (* d (conjugate c)))) eps)
      (incf flag 1000))

    (let* ((delta (+ (expt (- a d) 2) (* 4.0d0 b c)))
	   (deltaroot* (sqrt delta))
	   (deltaroot (if (minusp (realpart deltaroot*))
			  (- deltaroot*) deltaroot*))
	   (lambda1 (* 0.5d0 (+ a d deltaroot)))
	   (lambda2 (* 0.5d0 (+ a d (- deltaroot))))
	   (theta1 (- (argument lambda1)))
	   (theta2 (- (argument lambda2)))
	   ;; (theta0 (* 0.5d0 (+ theta1 theta2)))
	   beta
	   gamma
	   tmp)
      (when (> theta1 0.0d0)
	(psetf delta lambda1
	       lambda1 lambda2
	       lambda2 delta)
	;; Swap eigenvalues
	(psetf theta1 theta2 theta2 theta1))

      (when (> (abs (- (abs lambda1) 1.0d0)) eps)
	(incf flag 20))
      (when (> (abs (- (abs lambda2) 1.0d0)) eps)
	(incf flag 200))
      (setf beta (- (argument (- lambda1 a))
		    (argument b)))
      (when (< beta 0.0d0)
	(incf beta (* 2 pi)))
      (when (>= beta (* 2 pi))
	(decf beta (* 2 pi)))
      (setf gamma (acos (/ (- (expt (abs b) 2) (expt (abs (- lambda1 a)) 2))
			   (+ (expt (abs b) 2) (expt (abs (- lambda1 a)) 2)))))

      (setf tmp (sqrt (+ (expt (abs b) 2) (expt (abs (- lambda1 a)) 2))))
      (setf (aref eigenvectors 0 0) (/ b tmp))
      (setf (aref eigenvectors 1 0) (/ (- lambda1 a) tmp))
      (setf tmp (sqrt (+ (expt (abs b) 2) (expt (abs (- lambda2 a)) 2))))
      (setf (aref eigenvectors 0 1) (/ b tmp))
      (setf (aref eigenvectors 1 1) (/ (- lambda2 a) tmp))

      (setf deltaroot (+ (* (- a lambda1) (aref eigenvectors 0 0))
			 (* b (aref eigenvectors 1 0))))
      (setf delta (+ (* c (aref eigenvectors 0 0))
		     (* (- d lambda1) (aref eigenvectors 1 0))))
      (setf tmp (+ (abs deltaroot) (abs delta)))
      (when (> tmp eps)
	(incf flag 40))

      (setf deltaroot (+ (* (- a lambda2) (aref eigenvectors 0 1))
			 (* b (aref eigenvectors 1 1))))
      (setf delta (+ (* c (aref eigenvectors 0 1))
		     (* (- d lambda2) (aref eigenvectors 1 1))))
      (setf tmp (+ (abs deltaroot) (abs delta)))
      (when (> tmp eps)
	(incf flag 400))
      ;; Further Check ...
      (values theta1 theta2 gamma beta eigenvectors flag))))

;; Misc

#+nil
(defun sum (arr)
  (declare (type (vector double-float 5) arr))
  (loop
     for el across arr
     sum el))

(defun intercept-parameters (&rest params)
  "This is a dummy function existing solely for the purpose of being
   traceable with an arbitrary list of parameters."
  (declare (ignorable params)))
