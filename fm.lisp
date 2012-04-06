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

(defun store-cdfloat (fn seq)
  (with-open-file (s fn
		     :direction :output
		    :element-type '(unsigned-byte 8)
		    :if-does-not-exist :create
		    :if-exists :supersede)
   (let* ((a (sb-alien:sap-alien (sb-sys:vector-sap seq)
				 (sb-alien:* sb-alien:unsigned-char)))
	  (n (length seq))
	  (o (make-array (* n 16) :element-type '(unsigned-byte 8))))
     (dotimes (i (* 16 n))
       (setf (aref o i) (deref a i)))
     (write-sequence o s))))

(time
 (progn
  (defparameter *rate* 2048000)
  (defparameter *n-complex* (floor (expt 2 18) #+nil 47448064 2))
  (defparameter *input*
    (let* ((n (* 2 *n-complex*))
	   (a (make-array n :element-type '(unsigned-byte 8)))
	   (np2 (next-power-of-two (floor n 2)))
	   (c (make-array np2
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
		(* (exp (complex 0d0 (* (/ (* np2 592750) *rate*) (/ (* 2d0 pi) *n-complex*)
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

(defparameter *kin-filt*
  (let* ((n (length *kin*))
	 (nh (floor n 2))
	 
	 (bw (floor (* n 128000) ;; lowpass +/- 128kHz
		    *rate*))
	 (n-small (next-power-of-two (* 2 bw)))
	 (bww (floor n-small 2))
	 (res (make-array n-small :element-type (array-element-type *kin*)
			  :initial-element (complex .001d0))))
    (format t "2*bw=~6,2f bww=~d~%" bw bww)
    (loop for i from (+ (- bww) nh)
       below (+ bww nh) and ii from 0
       do
	 (setf (aref res ii) (aref *kin* i)))
    res))

(time
 (let* ((nn (length *kin-filt*))
       (all (make-array nn
			:element-type 'double-float))
	(m (map-into all #'abs *kin-filt*))
       (m2 m;(map-into all #'log all)
	 )
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

(time
 (progn
   (defparameter *input-filt*
     (napa-fft:fft
      (fftshift *kin-filt*)))
   (store-cdfloat "/dev/shm/kin-filt.cdfload"
		  *input-filt*)
   nil))


;; s = A e^ip
;; ds/dt = d/dt A e^ip = A ip' e^ip
;; ds/dt /s = ip'
;; p' = Im[ds/dt /s]
(time
 (progn
   (defparameter *demod-heterodyn*
    (let* ((in *input-filt*)
	   (n (length in))
	  
	   (d (make-array n 
			  :element-type 'double-float)))
      (loop for i from 1 below n do
	   (let* ((s (aref in i))
		  (ds/dt (- s (aref in (1- i)))))
	     (declare (type (complex double-float) s ds/dt))
	     (setf (aref d i) (imagpart (/ ds/dt
					   s)))))
      d))
       (store-dfloat "/dev/shm/demod-heterodyn.dfloat" *demod-heterodyn*)
       nil))

#+nil
(progn ;; print demodulated spectrum
 (let* ((k (napa-fft:fft *demod-heterodyn*))
	(mag (map-into (make-array (length k)
				   :element-type 'double-float)
		       #'abs (fftshift k))))
   (with-open-file (s "/dev/shm/o2.dat"
		      :direction :output
		      :if-does-not-exist :create
		      :if-exists :supersede)
     (loop for e across mag and i from (- (floor (length k) 2)) do
	  (format s "~6,3f ~6,3f~%"  (/ (* 2 i 128d3)
					(length k))  (log e))))))
(/ (* i 128000d0) (length k))
(* 3 1216)

(/ (* 2 (- (* 3 1216) 100) 128d3)
 (length *demod-heterodyn*))



(defparameter *rds*
 (progn ;; cut out +/-2kHz around 57kHz
   (let* ((bw 4000)
	  
	  (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	  (nn (length d))
	  (center-bin (floor (* 57d3 nn
				(/ 256d3))))
	  (n (next-power-of-two
	      (* 4000 nn
		 (/ (* 2 128d3)))))
	  (nh (floor n 2))
	  (a (make-array n :element-type (array-element-type d))))
     (loop for i from (- center-bin nh) below (+ center-bin nh)
	and ii from 0 
	do
	(setf (aref a ii) (aref d (+ (floor nn 2) i))))
     a)))

(defparameter *pilot*
 (progn ;; cut out +/-2kHz around 19kHz
   (let* ((bw 4000)
	  
	  (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	  (nn (length d))
	  (center-bin (floor (* 19d3 nn
				(/ 256d3))))
	  (n (next-power-of-two
	      (* 4000 nn
		 (/ (* 2 128d3)))))
	  (nh (floor n 2))
	  (a (make-array n :element-type (array-element-type d))))
     (loop for i from (- center-bin nh) below (+ center-bin nh)
	and ii from 0 
	do
	(setf (aref a ii) (aref d (+ (floor nn 2) i))))
     a)))

#+nil
(progn ;; print spectrum
 (let ((d *rds*))
   (with-open-file (s "/xdev/shm/o3.dat"
		      :direction :output
		      :if-does-not-exist :create
		      :if-exists :supersede)
     (loop for e across d and i from (- (floor (length d) 2)) do
	  (format s "~6,3f ~6,3f~%"  i  (abs e))))))


(napa-fft:fft 
 (fftshift *psk*))

(time
 (progn
   (defparameter *demod*
    (let* ((n (length *input-filt*		     ))
	   (d *input-filt*)
	   (e (make-array n 
			  :element-type 'double-float))
	   (o 1))
      (loop for j from o below n do
	   (let* ((ii (realpart (aref d (- j o))))
		  (i (realpart (aref d j)))
		  (q (imagpart (aref d j)))
		  (qq (imagpart (aref d (- j o)))))
	     (declare (type double-float i ii q qq))
	     (setf (aref e j) (- (* ii q) 
				 (* qq i)))))
      d))
   (store-dfloat "/dev/shm/demod.dfloat" *demod*)
   nil))

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
 )

