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
  (assert (every #'(lambda (x) (typep x type)) vectors))
  (assert (apply #'= (mapcar #'length vectors)))
  (macrolet ((vec-lin-com (coefficients vectors &key (type 'vector))
               (let ((args (loop for v in vectors collect (gensym))))
                 `(map ',type
                       (lambda ,args
                         (+ ,@(mapcar #'(lambda (c a) (list '* c a)) coefficients args)))
                       ,@vectors))))
    (vec-lin-com coefficients vectors :type type)))

;; (defun vector-linear-combination (coefficients vectors &key (type 'vector))
;;   (vec-lin-com coefficients vectors :type type))

;; (defmacro vec-lin-com (coefficients vectors &key (type 'vector))
;;   (let ((args (loop for v in vectors collect (gensym))))
;;     `(map ',type
;;           (lambda ,args
;;             (+ ,@(mapcar #'(lambda (c a) (list '* c a)) coefficients args)))
;;           ,@vectors)))

;; (vec-lin-com '(10 100 1000) '(#(1 2 3) #(4 5 6) #(7 8 9)))
