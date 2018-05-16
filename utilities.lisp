;;;; utilities.lisp

(in-package #:rf-utilities)

(defun sign (a b)
  "Implements a Fortran-like SIGN function."
  (* (abs a) (signum b)))

;; Matrix functions

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
  (loop
     for col from 0 below (array-dimension matrix 1)
     collect (aref matrix row col) into res
     finally (return (coerce res 'vector))))

(defun matrix-column (matrix col)
  (loop
     for row from 0 below (array-dimension matrix 0)
     collect (aref matrix row col) into res
     finally (return (coerce res 'vector))))

;; Complex functions

(setf (fdefinition 'argument) #'phase)
