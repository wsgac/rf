;;;; odedr-wanner-hairer.lisp

(in-package #:rf-odedr-wanner-hairer)

;;;; Methods for solving nonstiff systems of equations. Methods used
;;;; herein are described in Ernst Hairer, Gerhard Wanner and Syvert
;;;; P. Nørsett - Solving Ordinary Differential Equations I

;;;;;;;;;;;;;
;; General ;;
;;;;;;;;;;;;;

(defparameter uround double-float-epsilon)

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
  (let* ((ydopri1 (make-array n :initial-element 0))
         (kdopri1 (make-array n :initial-element 0))
         (kdopri2 (make-array n :initial-element 0))
         (kdopri3 (make-array n :initial-element 0))
         (kdopri4 (make-array n :initial-element 0))
         (kdopri5 (make-array n :initial-element 0))
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
    (loop
       do
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
         (setf h hnew))))

;;;;;;;;;;;;
;; DOPRI8 ;;
;;;;;;;;;;;;

;; DOPRI8 Coefficients
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
   H - guessed step of the independent variable"
  )

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
