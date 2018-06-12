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

;; Complex functions

(setf (fdefinition 'argument) #'phase)

;; Unitary operator functions

(defun angles-unitary (a b c d eps)
  (let ((flag 0)
        ())))

;; Misc

(defun sum (arr)
  (declare (type (vector double-float 5) arr))
  (loop
     for el across arr
     sum el))
