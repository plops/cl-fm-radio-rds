
(time
 (defparameter *input*
   (let* ((n 63438848)
	  (a (make-array n :element-type '(unsigned-byte 8)))
	  (c (make-array (floor n 2) :element-type '(complex single-float))))
     (with-open-file (s "/home/martin/Downloads2/dlf.fm"
			:element-type '(unsigned-byte 8))
       (read-sequence a s))
     (dotimes (i (floor n 2))
       (let ((x  (aref a (* 2 i)))
	     (y  (aref a (1+ (* 2 i)))))
	 (when (= x y)
	   (incf x .5))
	 (setf (aref c i) 
	       (complex (- x 127f0)
			(- y 127f0)))))
     c)))

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

#+nil
(let* ((n  (* 10 1024)) ;; verify data export with a sine wave
       (a (make-array n :element-type 'single-float)))
  (dotimes (i n)
    (setf (aref a i) (sin (* 2 #.(coerce pi 'single-float) i 200 (/ 1f0 n)))))
  (store-sfloat "/dev/shm/sin.sfloat" a)
  nil)

(time
 (store-sfloat "/dev/shm/demod.sfload" *demod*))