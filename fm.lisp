(ql:quickload "napa-fft3")

(defun next-power-of-two (n)
  (let ((p 2))
   (loop while (< p n) do 
	(setf p (* 2 p)))
   p))
#+nil
(next-power-of-two 33)

(defun store-sfloat (fn seq)
  (with-open-file (s fn
		     :direction :output
		    :element-type '(unsigned-byte 8)
		    :if-does-not-exist :create
		    :if-exists :supersede)
   (let* ((a (sb-alien:sap-alien (sb-sys:vector-sap seq)
				 (sb-alien:* sb-alien:unsigned-char)))
	  (n (length seq))
	  (o (make-array (* n 4) :element-type '(unsigned-byte 8))))
     (dotimes (i (* 4 n))
       (setf (aref o i) (deref a i)))
     (write-sequence o s))))

(defun store-dfloat (fn seq)
  (with-open-file (s fn
		     :direction :output
		    :element-type '(unsigned-byte 8)
		    :if-does-not-exist :create
		    :if-exists :supersede)
   (let* ((a (sb-alien:sap-alien (sb-sys:vector-sap seq)
				 (sb-alien:* sb-alien:unsigned-char)))
	  (n (length seq))
	  (o (make-array (* n 8) :element-type '(unsigned-byte 8))))
     (dotimes (i (* 8 n))
       (setf (aref o i) (deref a i)))
     (write-sequence o s))))

(time
 (progn
  (defparameter *rate* 2048000)
  (defparameter *n-complex* (floor (expt 2 16) #+nil 47448064 2))
  (defparameter *input*
    (let* ((n (* 2 *n-complex*))
	   (a (make-array n :element-type '(unsigned-byte 8)))
	   (c (make-array (next-power-of-two (floor n 2))
			  :element-type '(complex double-float))))
      (with-open-file (s "/home/martin/Downloads2/dl2048000.fm"
			 :element-type '(unsigned-byte 8))
	(read-sequence a s))
      (dotimes (i (floor n 2))
	(let ((x  (aref a (* 2 i)))
	      (y  (aref a (1+ (* 2 i)))))
	  (when (= x y)
	    (incf x .5))
	  (setf (aref c i) 
		(* #+nil (exp (complex 0d0 (* (/ (* 2d0 pi) *n-complex*)
					i)))
		   (complex (- x 127d0)
			    (- y 127d0))))))
      c))))


(defun fftshift (a)
  (let* ((n (length a))
	 (nh (floor n 2))
	 (b (make-array n
		       :element-type (array-element-type a)
		       :initial-element (complex .001d0))))
    (unless (evenp n)
      (error "n must be even."))
    (dotimes (i nh)
      (setf (aref b i) (aref a (- nh 1 i))
	    (aref b (+ nh i)) (aref a (- n 1 i))))
    b))


(time
 (defparameter *kin*
   (fftshift
    (napa-fft:fft *input*))))


(time
 (let* ((nn (length *kin*))
       (all (make-array nn
			:element-type 'double-float))
       (m (map-into all #'abs *kin*))
       (m2 (map-into all #'log all))
       (n 1000)
       (big-bins (floor nn
			n))
       (res (make-array n :element-type 'double-float
			:initial-element 0d0)))
  (defparameter *abs-kin* all
#+nil
    (dotimes (i n)
      (dotimes (b big-bins)
	(let ((p (+ b (* 1000 i))))
	  (when (< p nn)
	    (incf (aref res i) (aref all p)))))))
  (with-open-file (s "/dev/shm/o.dat"
		    :direction :output
		    :if-does-not-exist :create
		    :if-exists :supersede)
   (loop for i from (- (floor nn 2)) and e across *abs-kin* do
	(format s "~7,3f ~7,3f~%" (/ (* *rate* i) *n-complex*) e)))))

;; s = A e^ip
;; ds/dt = d/dt A e^ip = A ip' e^ip
;; ds/dt /s = ip'
;; p' = Im[ds/dt /s]
(time
 (defparameter *demod*
   (let* ((n (length *input*))
	  (d (make-array n 
			 :element-type 'single-float)))
     (loop for i from 1 below n do
	  (let* ((s (aref *input* i))
		 (ds/dt (- s (aref *input* (1- i)))))
	    (declare (type (complex single-float) s ds/dt))
	    (setf (aref d i) (imagpart (/ ds/dt
					  s)))))
     d)))

(time
 (defparameter *demod*
   (let* ((n (length *input*))
	  (d *input*)
	  (e (make-array n 
			 :element-type 'single-float))
	  (o 1))
     (loop for j from o below n do
	  (let* ((ii (realpart (aref d (- j o))))
		 (i (realpart (aref d j)))
		 (q (imagpart (aref d j)))
		 (qq (imagpart (aref d (- j o)))))
	    (declare (type single-float i ii q qq))
	    (setf (aref e j) (- (* ii q) 
				(* qq i)))))
     d)))

#+nil
(sb-ext:gc :full t)


#+nil
(let* ((n  (* 10 1024)) ;; verify data export with a sine wave
       (a (make-array n :element-type 'single-float)))
  (dotimes (i n)
    (setf (aref a i) (sin (* 2 #.(coerce pi 'single-float) i 200 (/ 1f0 n)))))
  (store-sfloat "/dev/shm/sin.sfloat" a)
  nil)

(time
 (store-sfloat "/dev/shm/demod.sfload" *demod*))

