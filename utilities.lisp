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
	(eigenvectors (make-array `(,rf::*n-equations* ,rf::*n-equations*)
				  :initial-element 0)))
    ;;     flag=0
    ;; tyci=abs(a)**2+abs(c)**2
    ;; if(abs(tyci-1.0_prec) > epsil)flag=flag+10
    ;; tyci=abs(b)**2+abs(d)**2
    ;; if(abs(tyci-1.0_prec) > epsil)flag=100+flag
    ;; tyci=abs(conjg(a)*b+conjg(c)*d)
    ;; if(tyci > epsil)flag=1000+flag
    (when (> (abs (- (+ (expt (abs a) 2) (expt (abs c) 2)) 1.0)) eps)
      (incf flag 10))
    (when (> (abs (- (+ (expt (abs b) 2) (expt (abs d) 2)) 1.0)) eps)
      (incf flag 100))
    (when (> (abs (+ (* b (conjugate a)) (* d (conjugate c)))) eps)
      (incf flag 1000))

;;         delta=(a-d)**2+4.0_prec*b*c
;;     deltaroot=sqrt(delta)
;;     if(real(deltaroot,prec) < 0.0_prec)deltaroot=-deltaroot
;;     lambda1=0.5_prec*(a+d+deltaroot)
;;     lambda2=0.5_prec*(a+d-deltaroot)
;; !!$    wyznaczam fazy theta1,theta2
;;     theta1=-fargument(lambda1)
;;     theta2=-fargument(lambda2)
;;     if(theta1 > 0.0_prec)then !zamieniam miejscami wartosci wlasne
;;        delta=lambda1
;;        lambda1=lambda2
;;        lambda2=delta
;;        tyci=theta1
;;        theta1=theta2
;;        theta2=tyci
;;     endif
;;     theta0=0.5_prec*(theta1+theta2)
    (let* ((delta (+ (expt (- a d) 2) (* 4.0d0 b c)))
	   (deltaroot* (sqrt delta))
	   (deltaroot (if (minusp (realpart deltaroot*))
			  (- deltaroot*) deltaroot*))
	   (lambda1 (* 0.5d0 (+ a d deltaroot)))
	   (lambda2 (* 0.5d0 (+ a d (- deltaroot))))
	   (theta1 (- (argument lambda1)))
	   (theta2 (- (argument lambda2)))
	   (theta0 (* 0.5d0 (+ theta1 theta2)))
	   beta
	   gamma
	   tmp)
      (when (> theta1 0.0d0)
	(psetf delta lambda1
	       lambda1 lambda2
	       lambda2 delta)
	;; Swap eigenvalues
	(psetf theta1 theta2 theta2 theta1))

;;           tyci=abs(lambda1)
;;     if(abs(tyci-1.0_prec) > epsil)flag=flag+20
;;     tyci=abs(lambda2)
;;     if(abs(tyci-1.0_prec) > epsil)flag=200+flag
;; !!$    if(flag > 0)then
;; !!$       print*,'wartosci wlasne nie sa czynnikami fazowymi, flag=',flag
;; !!$       stop
;; !!$    endif
;; !!$    ===================================
;; !!$    flag=0
;; !!$    wyznaczam faze beta
;;     beta=fargument(lambda1-a)-fargument(b)
;;     if(beta < 0.0_prec) beta=beta+dwapi
;;     if(beta >= dwapi) beta=beta-dwapi
;; !!$    wyznaczam faze gama
;;     gama=(abs(b)**2-abs(lambda1-a)**2)/(abs(b)**2+abs(lambda1-a)**2)
;;     gama=acos(gama)
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

;;       !!$    z dokladnoscia do globalnej fazy
;;     tyci=sqrt(abs(b)**2+abs(lambda1-a)**2)
;;     eigenwektory(1,1)= b/tyci
;;     eigenwektory(2,1)= (lambda1-a)/tyci
;;     tyci=sqrt(abs(b)**2+abs(lambda2-a)**2)
;;     eigenwektory(1,2)= b/tyci
;;     eigenwektory(2,2)= (lambda2-a)/tyci
;; !!$    sprawdzam warunek wlasny
;;     deltaroot=(a-lambda1)*eigenwektory(1,1)+b*eigenwektory(2,1)
;;     delta    =c*eigenwektory(1,1)+(d-lambda1)*eigenwektory(2,1)
;;     tyci=abs(deltaroot)+abs(delta)
;;     if(tyci > epsil)flag=40+flag
;;     deltaroot=(a-lambda2)*eigenwektory(1,2)+b*eigenwektory(2,2)
;;     delta    =c*eigenwektory(1,2)+(d-lambda2)*eigenwektory(2,2)
;;     tyci=abs(deltaroot)+abs(delta)
;;     if(tyci > epsil)flag=400+flag
;;     iGlobalCheck=flag
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
	(incf flag 400)))))

;; Misc

(defun sum (arr)
  (declare (type (vector double-float 5) arr))
  (loop
     for el across arr
     sum el))
