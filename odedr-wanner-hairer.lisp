;;;; odedr-wanner-hairer.lisp

(in-package #:rf-odedr-wanner-hairer)

;;;; Methods for solving nonstiff systems of equations. Methods used
;;;; herein are described in Ernst Hairer, Gerhard Wanner and Syvert
;;;; P. Nørsett - Solving Ordinary Differential Equations I

;;;;;;;;;;;;;
;; General ;;
;;;;;;;;;;;;;

(defparameter uround double-float-epsilon)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; ;;;;;;;;;;;;;;;;;;;;; ;;
;; ;; Solving Methods ;; ;;
;; ;;;;;;;;;;;;;;;;;;;;; ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;
;; DOPRI5 ;;
;;;;;;;;;;;;

;; DOPRI5 Coefficients
(defparameter nmax_dopri5 10000)
(defparameter c2 (/ 1.0 18.0))          ; Taken from DOPRI8
(defparameter cc2 (/ 1.0 5.0))
(defparameter cc3 (/ 3.0 10.0))
(defparameter cc4 (/ 4.0 5.0))
(defparameter cc5 (/ 8.0 9.0))
(defparameter aa21 c2)
(defparameter aa31 (/ 3.0 40.0))
(defparameter aa32 (/ 9.0 40.0))
(defparameter aa41 (/ 44.0 45.0))
(defparameter aa42 (/ -56.0 15.0))
(defparameter aa43 (/ 32.0 9.0))
(defparameter aa51 (/ 19372.0 6561.0))
(defparameter aa52 (/ -25360.0 2187.0))
(defparameter aa53 (/ 64448.0 6561.0))
(defparameter aa54 (/ -212.0 729.0))
(defparameter aa61 (/ 9017.0 3168.0))
(defparameter aa62 (/ -355.0 33.0))
(defparameter aa63 (/ 46732.0 5247.0))
(defparameter aa64 (/ 49.0 176.0))
(defparameter aa65 (/ -5103.0 18656.0))
(defparameter bb1 (/ 35.0 384.0))
(defparameter bb2 (/ 500.0 1113.0))
(defparameter bb3 (/ 125.0 192.0))
(defparameter bb4 (/ -2187.0 6784.0))
(defparameter bb5 (/ 11.0 84.0))
(defparameter bb6 (/ 71.0 57600.0))
(defparameter bb7 (/ -71.0 16695.0))
(defparameter bb8 (/ 71.0 1920.0))
(defparameter bb9 (/ -17253.0 339200.0))
(defparameter bb10 (/ 22.0 525.0))
(defparameter dd1 (/ -1.0 40.0))


(defun dopri5 (n fn x y xend epsilon hmax h ;iflagdopri5
               )
  "Implement the DOPRI5 method of solving systems of ODEs

   Parameter description:
   
   N - number of equations in a system
   FN - function of the system of equations (must take the arguments:
   N, X, Y, YP and return the modified YP vector)
   X - initial value of the independent variable
   Y - array (of size N) of initial values of solutions
   XEND - final value of the independent variable
   EPSILON - precision of the solutions
   HMAX - maximum step of the independent variable
   H - guessed step of the independent variable
   ;;IFLAGDOPRI5 - has XEND been reached after a specific number of steps
   (0 - done; 1 - not done, continue integration; 2 - probably a stiff
   system, change method)"
  (prog* (;; Auxiliary arrays
          (ydopri1 (make-array n :initial-element 0))
          (kdopri1 (make-array n :initial-element 0))
          (kdopri2 (make-array n :initial-element 0))
          (kdopri3 (make-array n :initial-element 0))
          (kdopri4 (make-array n :initial-element 0))
          (kdopri5 (make-array n :initial-element 0))
          ;; Remaining variables
          (iflagdopri5 0)
          (posneg (sign 1.0 (- xend x)))
          (hmax (abs hmax))
          (h (min (max 1.0e-4 (abs h)) hmax))
          (h (sign h posneg))
          (epsilon (max epsilon (* 7.0 uround)))
          (reject_dopri5 nil)
          (xph 0.0)
          (err 0.0)
          (denom 0.0)
          (fac 0.0)
          (hnew 0.0)
          (nstep 0)
          (nfcn 0)
          (naccept 0)
          (nreject 0))
   point1
   (when (> nstep nmax_dopri5)
     (setf iflagdopri5 1)
     (return (values iflagdopri5)))
   (when (= x (+ x (* 0.1 h)))
     (setf iflagdopri5 2)
     (return (values iflagdopri5)))
   (when (> (+ uround (* posneg (- x xend))) 0.0)
     (return (values iflagdopri5)))
   (when (> (* posneg (+ x h (- xend))) 0.0)
     (setf h (- xend x)))
   ;; Series of FN invocations
   (setf kdopri1 (funcall fn n x y kdopri1))
   (incf nstep)
   (setf ydopri1 (map 'vector
                      (lambda (el) (+ y (* h aa21 el)))
                      kdopri1))
   (setf kdopri2 (funcall fn n (+ x (* cc2 h)) ydopri1 kdopri2))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2) (+ y (* h (+ (* aa31 el1) (* aa32 el2)))))
                      kdopri1 kdopri2))
   (setf kdopri3 (funcall fn n (+ x (* c3 h)) ydopri1 kdopri3))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3)
                        (+ y (* h (+ (* aa41 el1) (* aa42 el2) (* aa43 el3)))))
                      kdopri1 kdopri2 kdopri3))
   (setf kdopri4 (funcall fn n (+ x (* c4 h)) ydopri1 kdopri4))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4)
                        (+ y (* h (+ (* aa51 el1) (* aa52 el2) (* aa53 el3) (* aa54 el4)))))
                      kdopri1 kdopri2 kdopri3 kdopri4))
   (setf kdopri5 (funcall fn n (+ x (* c5 h)) ydopri1 kdopri5))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4 el5)
                        (+ y (* h (+ (* aa61 el1) (* aa62 el2) (* aa63 el3) (* aa64 el4) (* aa65 el5)))))
                      kdopri1 kdopri2 kdopri3 kdopri4 kdopri5))
   (setf xph (+ x h))
   (setf kdopri2 (funcall fn n xph ydopri1 kdopri2))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4 el5)
                        (+ y (* h (+ (* bb1 el1) (* bb2 el2) (* bb3 el3) (* bb4 el4) (* bb5 el5)))))
                      kdopri1 kdopri3 kdopri4 kdopri5 kdopri2))
   (setf kdopri2 (map 'vector
                      (lambda (el1 el2 el3 el4 el5)
                        (+ (* bb6 el1) (* bb7 el2) (* bb8 el3) (* bb9 el4) (* bb10 el5)))
                      kdopri1 kdopri3 kdopri4 kdopri5 kdopri2))
   (setf kdopri3 (funcall fn n xph ydopri1 kdopri3))
   (setf kdopri4 (map 'vector
                      (lambda (el1 el2)
                        (* h (+ el1 (* dd1 el2))))
                      kdopri2 kdopri3))
   (incf nfcn 6)
   (setf err 0.0)
   (loop
      for i from 0 below n
      do (setf denom (max 1.0e-5
                          (abs (svref ydopri1 i))
                          (abs (svref y i))
                          (/ (* 0.2 uround) epsilon)))
        (incf err (expt (/ (svref kdopri4 i) denom) 2)))
   (setf err (sqrt (/ err n)))
   (setf fac (max 0.1
                  (min 5.0
                       (/ (expt (/ err epsilon) 0.2) 0.9))))
   (setf hnew (/ h fac))
   ;; Final iteration check
   (if (< err epsilon)
       (progn
         (incf naccept)
         (setf kdopri1 kdopri3)
         (setf y ydopri1)
         (setf x xph)
         (when (> (abs hnew) hmax)
           (setf hnew (* posneg hmax)))
         (when reject_dopri5
           (setf hnew (* posneg (min (abs hnew) (abs h)))))
         (setf reject_dopri5 nil))
       (progn
         (setf reject_dopri5 t)
         (when (> naccept 1)
           (incf nreject))))
   (setf h hnew)
   (go point1)))

;;;;;;;;;;;;
;; DOPRI8 ;;
;;;;;;;;;;;;

;; DOPRI8 Coefficients
(defparameter nmax_dopri8 100000000)
;; (defparameter c2 (/ 1.0 18.0))
(defparameter c3 (/ 1.0 12.0))
(defparameter c4 (/ 1.0 8.0))
(defparameter c5 (/ 5.0 16.0))
(defparameter c6 (/ 3.0 8.0))
(defparameter c7 (/ 59.0 400.0))
(defparameter c8 (/ 93.0 200.0))
(defparameter c9 (/ 5490023248.0 9719169821.0))
(defparameter c10 (/ 13.0 20.0))
(defparameter c11 (/ 1201146811.0 1299019798.0))
(defparameter c12 1.0)
(defparameter c13 1.0)
(defparameter a21 c2)
(defparameter a31 (/ 1.0 48.0))
(defparameter a32 (/ 1.0 16.0))
(defparameter a41 (/ 1.0 32.0))
(defparameter a43 (/ 3.0 32.0))
(defparameter a51 (/ 5.0 16.0))
(defparameter a53 (/ -75.0 64.0))
(defparameter a54 (- a53))
(defparameter a61 (/ 3.0 80.0))
(defparameter a64 (/ 3.0 16.0))
(defparameter a65 (/ 3.0 20.0))
(defparameter a71 (/ 29443841.0 614563906.0))
(defparameter a74 (/ 77736538.0 692538347.0))
(defparameter a75 (/ -28693883.0 1125.0e6))
(defparameter a76 (/ 23124283.0 18.0e8))
(defparameter a81 (/ 16016141.0 946692911.0))
(defparameter a84 (/ 61564180.0 158732637.0))
(defparameter a85 (/ 22789713.0 633445777.0))
(defparameter a86 (/ 545815736.0 2771057229.0))
(defparameter a87 (/ -180193667.0 1043307555.0))
(defparameter a91 (/ 39632708.0 573591083.0))
(defparameter a94 (/ -433636366.0 683701615.0))
(defparameter a95 (/ -421739975.0 2616292301.0))
(defparameter a96 (/ 100302831.0 723423059.0))
(defparameter a97 (/ 790204164.0 839813087.0))
(defparameter a98 (/ 800635310.0 3783071287.0))
(defparameter a101 (/ 246121993.0 1340847787.0))
(defparameter a104 (/ -37695042795.0 15268766246.0))
(defparameter a105 (/ -309121744.0 1061227803.0))
(defparameter a106 (/ -12992083.0 490766935.0))
(defparameter a107 (/ 6005943493.0 2108947869.0))
(defparameter a108 (/ 393006217.0 1396673457.0))
(defparameter a109 (/ 123872331.0 1001029789.0))
(defparameter a111 (/ -1028468189.0 846180014.0))
(defparameter a114 (/ 8478235783.0 508512852.0))
(defparameter a115 (/ 1311729495.0 1432422823.0))
(defparameter a116 (/ -10304129995.0 1701304382.0))
(defparameter a117 (/ -48777925059.0 3047939560.0))
(defparameter a118 (/ 15336726248.0 1032824649.0))
(defparameter a119 (/ -45442868181.0 3398467696.0))
(defparameter a1110 (/ 3065993473.0 597172653.0))
(defparameter a121 (/ 185892177.0 718116043.0))
(defparameter a124 (/ -3185094517.0 667107341.0))
(defparameter a125 (/ -477755414.0 1098053517.0))
(defparameter a126 (/ -703635378.0 230739211.0))
(defparameter a127 (/ 5731566787.0 1027545527.0))
(defparameter a128 (/ 5232866602.0 850066563.0))
(defparameter a129 (/ -4093664535.0 808688257.0))
(defparameter a1210 (/ 3962137247.0 1805957418.0))
(defparameter a1211 (/ 65686358.0 487910083.0))
(defparameter a131 (/ 403863854.0 491063109.0))
(defparameter a134 (/ -5068492393.0 434740067.0))
(defparameter a135 (/ -411421997.0 543043805.0))
(defparameter a136 (/ 652783627.0 914296604.0))
(defparameter a137 (/ 11173962825.0 925320556.0))
(defparameter a138 (/ -13158990841.0 6184727034.0))
(defparameter a139 (/ 3936647629.0 1978049680.0))
(defparameter a1310 (/ -160528059.0 685178525.0))
(defparameter a1311 (/ 248638103.0 1413531060.0))
(defparameter b1 (/ 14005451.0 335480064.0))
(defparameter b6 (/ -59238493.0 1068277825.0))
(defparameter b7 (/ 181606767.0 758867731.0))
(defparameter b8 (/ 561292985.0 797845732.0))
(defparameter b9 (/ -1041891430.0 1371343529.0))
(defparameter b10 (/ 760417239.0 1151165299.0))
(defparameter b11 (/ 118820643.0 751138087.0))
(defparameter b12 (/ -528747749.0 2220607170.0))
(defparameter b13 (/ 1.0 4.0))
(defparameter bh1 (/ 13451932.0 455176623.0))
(defparameter bh6 (/ -808719846.0 976000145.0))
(defparameter bh7 (/ 1757004468.0 5645159321.0))
(defparameter bh8 (/ 656045339.0 265891186.0))
(defparameter bh9 (/ -3867574721.0 1518517206.0))
(defparameter bh10 (/ 465885868.0 322736535.0))
(defparameter bh11 (/ 53011238.0 667516719.0))
(defparameter bh12 (/ 2.0 45.0))

(defun dopri8 (n fn x y xend epsilon hmax h)
  "Implement the DOPRI8 method of solving systems of ODEs

   Parameter description:
   
   N - number of equations in a system
   FN - function of the system of equations (must take the arguments:
   N, X, Y, YP and return the modified YP vector)
   X - initial value of the independent variable
   Y - array (of size N) of initial values of solutions
   XEND - final value of the independent variable
   EPSILON - precision of the solutions
   HMAX - maximum step of the independent variable
   H - guessed step of the independent variable

   This is the implementation of a real-number variant of DOPRI8."
  (prog* ( ;; Auxiliary arrays
          (ydopri1 (make-array n :initial-element 0))
          (kdopri1 (make-array n :initial-element 0))
          (kdopri2 (make-array n :initial-element 0))
          (kdopri3 (make-array n :initial-element 0))
          (kdopri4 (make-array n :initial-element 0))
          (kdopri5 (make-array n :initial-element 0))
          (kdopri6 (make-array n :initial-element 0))
          (kdopri7 (make-array n :initial-element 0))
          (y11s (make-array n :initial-element 0))
          (y12s (make-array n :initial-element 0))
          ;; Remaining variables
          (iflagdopri8 0)
          (posneg (sign 1.0 (- xend x)))
          (hmax (abs hmax))
          (h (min (max 1.0e-10 (abs h)) hmax))
          (h (sign h posneg))
          (epsilon (max epsilon (* 13.0 uround)))
          (reject_dopri8 nil)
          (xph 0.0)
          (err 0.0)
          (denom 0.0)
          (fac 0.0)
          (hnew 0.0)
          (nstep 0)
          (nfcn 0)
          (naccept 0)
          (nreject 0))

   point1 ;; Equivalent to 1

   (when (zerop iflagdopri8)
     (when (> nstep nmax_dopri8)
       (setf iflagdopri8 1)
       (return (values y iflagdopri8)))
     (when (= x (+ x (* 0.03 h)))
       (setf iflagdopri8 2)
       (return (values y iflagdopri8))))
   (when (> (+ uround (* posneg (- x xend))) 0.0)
     (return (values y iflagdopri8)))
   (when (> (* posneg (+ x h (- xend))) 0.0)
     (setf h (- xend x)))
   (setf kdopri1 (funcall fn n x y kdopri1))

   point2 ;; Equivalent to 2

   (incf nstep)
   ;; Sequence of function calls
   (setf ydopri1 (map 'vector
                      (lambda (el1)
                        (+ y (* h a21 el1)))
                      kdopri1))
   (setf kdopri2 (funcall fn n (+ x (* c2 h)) ydopri1 kdopri2))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2)
                        (+ y (* h (+ (* a31 el1) (* a32 el2)))))
                      kdopri1 kdopri2))
   (setf kdopri3 (funcall fn n (+ x (* c3 h)) ydopri1 kdopri3))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2)
                        (+ y (* h (+ (* a41 el1) (* a43 el2)))))
                      kdopri1 kdopri3))
   (setf kdopri4 (funcall fn n (+ x (* c4 h)) ydopri1 kdopri4))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3)
                        (+ y (* h (+ (* a51 el1) (* a53 el2) (* a54 el3)))))
                      kdopri1 kdopri3 kdopri4))
   (setf kdopri5 (funcall fn n (+ x (* c5 h)) ydopri1 kdopri5))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3)
                        (+ y (* h (+ (* a61 el1) (* a64 el2) (* a65 el3)))))
                      kdopri1 kdopri4 kdopri5))
   (setf kdopri6 (funcall fn n (+ x (* c6 h)) ydopri1 kdopri6))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4)
                        (+ y (* h (+ (* a71 el1) (* a74 el2) (* a75 el3) (* a76 el4)))))
                      kdopri1 kdopri4 kdopri5 kdopri6))
   (setf kdopri7 (funcall fn n (+ x (* c7 h)) ydopri1 kdopri7))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4 el5)
                        (+ y (* h (+ (* a81 el1) (* a84 el2) (* a85 el3) (* a86 el4) (* a87 el5)))))
                      kdopri1 kdopri4 kdopri5 kdopri6 kdopri7))
   (setf kdopri2 (funcall fn n (+ x (* c8 h)) ydopri1 kdopri2))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4 el5 el6)
                        (+ y (* h (+ (* a91 el1) (* a94 el2) (* a95 el3) (* a96 el4) (* a97 el5) (* a98 el6)))))
                      kdopri1 kdopri4 kdopri5 kdopri6 kdopri7 kdopri2))
   (setf kdopri3 (funcall fn n (+ x (* c9 h)) ydopri1 kdopri3))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3 el4 el5 el6 el7)
                        (+ y (* h (+ (* a101 el1) (* a104 el2) (* a105 el3) (* a106 el4) (* a107 el5) (* a108 el6) (* a109 el7)))))
                      kdopri1 kdopri4 kdopri5 kdopri6 kdopri7 kdopri2 kdopri3))
   ;; Some variables
   ;; y11s=a111*kdopri1+a114*kdopri4+a115*kdopri5+a116*kdopri6+a117*kdopri7+a118*kdopri2+a119*kdopri3
   ;; y12s=a121*kdopri1+a124*kdopri4+a125*kdopri5+a126*kdopri6+a127*kdopri7+a128*kdopri2+a129*kdopri3
   ;; kdopri4=a131*kdopri1+a134*kdopri4+a135*kdopri5+a136*kdopri6+a137*kdopri7+a138*kdopri2+a139*kdopri3
   ;; kdopri5=b1*kdopri1+b6*kdopri6+b7*kdopri7+b8*kdopri2+b9*kdopri3
   ;; kdopri6=bh1*kdopri1+bh6*kdopri6+bh7*kdopri7+bh8*kdopri2+bh9*kdopri3
   ;; kdopri2=y11s
   ;; kdopri3=y12s
   (setf y11s (vector-linear-combination (list a111 a114 a115 a116 a117 a118 a119)
                                         (list kdopri1 kdopri4 kdopri5 kdopri6 kdopri7 kdopri2 kdopri3)))
   (setf y12s (vector-linear-combination (list a121 a124 a125 a126 a127 a128 a129)
                                         (list kdopri1 kdopri4 kdopri5 kdopri6 kdopri7 kdopri2 kdopri3)))
   (setf kdopri4 (vector-linear-combination (list a131 a134 a135 a136 a137 a138 a139)
                                            (list kdopri1 kdopri4 kdopri5 kdopri6 kdopri7 kdopri2 kdopri3)))
   (setf kdopri5 (vector-linear-combination (list b1 b6 b7 b8 b9)
                                            (list kdopri1 kdopri6 kdopri7 kdopri2 kdopri3)))
   (setf kdopri6 (vector-linear-combination (list bh1 bh6 bh7 bh8 bh9)
                                            (list kdopri1 kdopri6 kdopri7 kdopri2 kdopri3)))
   (setf kdopri2 y11s)
   (setf kdopri3 y12s)
   ;; Another sequence of function calls
   (setf kdopri7 (funcall fn n (+ x (* c10 h)) ydopri1 kdopri7))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2)
                        (+ y (* h (+ (* 1 el1) (* a1110 el2)))))
                      kdopri2 kdopri7))
   (setf kdopri2 (funcall fn n (+ x (* c11 h)) ydopri1 kdopri2))
   (setf xph (+ x h))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3)
                        (+ y (* h (+ (* 1 el1) (* a1210 el2) (* a1211 el3)))))
                      kdopri3 kdopri7 kdopri2))
   (setf kdopri3 (funcall fn n xph ydopri1 kdopri3))
   (setf ydopri1 (map 'vector
                      (lambda (el1 el2 el3)
                        (+ y (* h (+ (* 1 el1) (* a1310 el2) (* a1311 el3)))))
                      kdopri4 kdopri7 kdopri2))
   (setf kdopri4 (funcall fn n xph ydopri1 kdopri4))
   (incf nfcn 13)
   (setf kdopri5 (vector-linear-combination (list 1 b10 b11 b12 b13)
                                            (list kdopri5 kdopri7 kdopri2 kdopri3 kdopri4)))
   (setf kdopri6 (vector-linear-combination (list 1 bh10 bh11 bh12)
                                            (list kdopri6 kdopri7 kdopri2 kdopri3)))

   (setf err 0.0)
   (loop
      for i from 0 below n
      do (setf denom (max 1.0e-6
                          (abs (svref kdopri5 i))
                          (abs (svref y i))
                          (/ (* 0.2 uround) epsilon)))
        (incf err (expt (/ (- (svref kdopri5 i) (svref kdopri6 i)) denom) 2)))
   (setf err (sqrt (/ err n)))

   (setf fac (max (/ 1.0 6.0) (min 3.0 (/ (expt (/ err epsilon) (/ 1.0 8.0)) 0.9))))
   (setf hnew (/ h fac))
   (when (> err epsilon)
     (go point51))                      ; GOTO
   (incf naccept)
   (setf y kdopri5)

   (setf x xph)
   (when (> (abs hnew) hmax)
     (setf hnew (* posneg hmax)))
   (when reject_dopri8
     (setf hnew (* posneg (min (abs hnew) (abs h)))))
   (setf reject_dopri8 nil)
   (setf h hnew)
   (go point1)                          ; GOTO
   
   point51 ;; Equivalent to 51

   (setf reject_dopri8 t)
   (setf h hnew)
   (when (> naccept 1)
     (incf nreject))
   (decf nfcn)
   (go point2)))                        ; GOTO

;;;;;;;;;;
;; ODEX ;;
;;;;;;;;;;

;; ODEX Coefficients
(defparameter a #(3.0 7.0 13.0 21.0 31.0 43.0 57.0 73.0 91.0))
(defparameter nj #(2 4 6 8 10 12 14 16 18))
(defparameter fac1 2.0e-2)
(defparameter fac2 4.0)
(defparameter fac3 0.9)
(defparameter fac4 0.8)
(defparameter safe1 0.65)
(defparameter safe2 0.94)

;;;;;;;;;;;;;;;;;;;;;;
;; Monodromy Matrix ;;
;;;;;;;;;;;;;;;;;;;;;;

(defun monodromy-matrix (xin xend epsilonabs e-amplitude sigma n-cycles k-perpendicular cep basic-shape)
  )
