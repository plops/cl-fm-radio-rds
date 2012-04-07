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
  (defparameter *n-complex* (floor (expt 2 25) #+nil 47448064 2))
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
		(* (exp (complex 0d0 (* (/ (* np2 -507250 ;592750
					      ) *rate*) (/ (* 2d0 pi) *n-complex*)
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
	 
	 (bw (floor (* n 128000) ;; lowpass +/- 64kHz
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
 (progn ;; print full spectrum
  (let* ((nn (length *kin-filt*))
	 (all (make-array nn
			  :element-type 'double-float))
	 (m (map-into all #'abs *kin-filt*))
	 (m2 m				;(map-into all #'log all)
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
	   (format s "~7,3f ~7,3f~%" (/ (* *rate* i) *n-complex*) e))))))

(time
 (progn ;; reverse ft
   (defparameter *input-filt*
     (napa-fft:ifft
      (fftshift *kin-filt*)))
   (store-cdfloat "/dev/shm/kin-filt.cdfload"
		  *input-filt*)
   nil))


;; s = A e^ip
;; ds/dt = d/dt A e^ip = A ip' e^ip
;; ds/dt /s = ip'
;; p' = Im[ds/dt /s]
(time
 (progn ;; demodulate using heterodyne division
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

(defparameter *pilot*
 (progn ;; cut out +/-250Hz around 19kHz
   (let* ((bw 500)
	  (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	  (nn (length d))
	  (center-bin (floor (* 19d3 nn
				(/ 256d3))))
	  (band-bin (* bw nn
		       (/ (* 2 128d3))))
	  (nh (floor band-bin 2))
	  (a (make-array nn :element-type (array-element-type d))))
     (loop for i from (- center-bin nh) below (+ center-bin nh)
	do
	  (setf (aref a (+ (floor nn 2) i)) (aref d (+ (floor nn 2) i))
		(aref a (- (floor nn 2) i)) (aref d (- (floor nn 2) i))))
     a)))

(defparameter *pilot3*
 (progn ;; cut out +/-1000Hz around 57kHz
   (let* ((bw 2000)
	  (rp (napa-fft:fft (fftshift *pilot*)))
	  (ma (reduce #'(lambda (x y) (max (realpart x) (realpart y))) (subseq rp 0 4000)))
	  (rp2 (map-into rp #'(lambda (z) (complex (expt (/ (realpart z) ma) 3))) rp))
	  (d (fftshift (napa-fft:fft rp2)))
	  (nn (length d))
	  (center-bin (floor (* 57d3 nn
				(/ 256d3))))
	  (band-bin (* bw nn
		       (/ (* 2 128d3))))
	  (nh (floor band-bin 2))
	  (a (make-array nn :element-type (array-element-type d))))
     (loop for i from (- center-bin nh) below (+ center-bin nh)
	do
	  (setf (aref a (+ (floor nn 2) i)) (aref d (+ (floor nn 2) i))
		(aref a (- (floor nn 2) i)) (aref d (- (floor nn 2) i))))
     a)))

(defparameter *rds*
 (progn ;; cut out +/-2kHz around 57kHz
   (let* ((bw 4000)
	  (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	  (nn (length d))
	  (center-bin (floor (* 57d3 nn
				(/ 256d3))))
	  (n (* bw nn
		(/ (* 2 128d3))))
	  (nh (floor n 2))
	  (a (make-array nn :element-type (array-element-type d))))
     (loop for i from (- center-bin nh) below (+ center-bin nh)
	do
	  (setf (aref a (+ (floor nn 2) i)) (aref d (+ (floor nn 2) i))
		(aref a (- (floor nn 2) i)) (aref d (- (floor nn 2) i))))
     a)))

(progn ;; 
 (defparameter *rds-base*
   (progn ;; cut out +/-2kHz around 57kHz
     (let* ((bw 4000)
	    (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	    (nn (length d))
	    (center-bin (floor (* 57d3 nn
				  (/ 256d3))))
	    (n (floor (* bw nn
			 (/ (* 2 128d3)))))
	    (nh (floor n 2))
	    (a (make-array n :element-type (array-element-type d))))
       (loop for i from (- center-bin nh) below (+ center-bin nh)
	  and ii from 0
	  do
	  (setf (aref a ii) (aref d (+ (floor nn 2) i))))
       a)))
 (store-cdfloat "/dev/shm/rds-base-square.cdfloat" (let* ((a (make-array (length *demod-heterodyn*)
								  :element-type '(complex double-float)))
						   (b (map-into a #'(lambda (x) (* x x))
								(napa-fft:ifft (fftshift *rds-base*)))))
					      b))
 nil)

(floor (* 57d3 (length *demod-heterodyn*)
	  (/ 256d3))
       48)
9728
2097152
32768

(* 32768 (/ 9728 2097152))
;; average 152 values


(length *pilot*)

(* 3 156000)

(/ 19000 48f0)

(/ 256000 1187.5)
#+nil
(progn ;; print spectrum
 (let ((d *pilot3*))
   (with-open-file (s "/dev/shm/o3.dat"
		      :direction :output
		      :if-does-not-exist :create
		      :if-exists :supersede)
     (loop for e across d and i from (- (floor (length d) 2)) do
	  (format s "~6,3f ~6,3f~%"  i  (log (abs (+ .0001 e))))))))

(progn ;; print waveform
 (let ((d (napa-fft:fft (fftshift *pilot*)))
       (d3 (napa-fft:fft (fftshift *pilot3*)))
       (f (napa-fft:fft (fftshift *rds*))))
   (with-open-file (s "/dev/shm/o3.dat"
		      :direction :output
		      :if-does-not-exist :create
		      :if-exists :supersede)
     #+nil
     (loop for pil across d and sig across f and i from 0 below 30000 do
	  (let ((sum (complex 0d0)))
	    (incf sum (complex (realpart pil) (realpart sig)))
	    (when (= 0 (mod (+ 215 i) 20))
	      (format s "~6,3f ~6,3f~%" (realpart sum) (imagpart sum))
	      (setf sum (complex 0d0)
		))))

     #+nil
     (progn ;; mul
       (terpri s)
       (loop for pil across d3 and sig across f and i from 0 below 4000 do
	    (let ((q (complex (realpart pil)
			      (realpart sig))))
	      (setf q (/ q (abs q)))
	      (format s "~6,3f ~6,3f~%"  i (let ((p (let ((v (realpart pil)))
						     v #+nil (cond ((< v -10) -200)
							    ((<= v 10) v)
							    ((< 10 v) 200))))
						 (q (let ((v (realpart sig)))
						      v #+nil(cond ((< v -10) -200)
							    ((<= v 10) v)
							    ((< 10 v) 200)))))
					     (* p q))
		     ))))
     (progn ;; pilot3 rect
       (terpri s)
       (loop for pil across d3 and sig across f and i from 0 below 40000 do
	    (let ((q (complex (realpart pil)
			      (realpart sig))))
	      (setf q (/ q (abs q)))
	      (format s "~6,3f ~6,3f~%"  (realpart sig) (let ((v (realpart pil)))
					     v #+nil (cond ((< v -10) -200000)
						   ((<= v 10) v)
						   ((< 10 v) 200000)))
		     ))))
     #+nil(progn ;; pilot rect
       (terpri s)
       (loop for pil across d and sig across f and i from 0 below 4000 do
	    (let ((q (complex (realpart pil)
			      (realpart sig))))
	      (setf q (/ q (abs q)))
	      (format s "~6,3f ~6,3f~%"  i (let ((v (realpart pil)))
					     (cond ((< v -10) -300000)
						   ((<= v 10) v)
						   ((< 10 v) 300000)))
		     ))))
     #+nil(progn ;; rds rect
       (terpri s)
       (loop for pil across d and sig across f and i from 0 below 4000 do
	    (let ((q (complex (realpart pil)
			      (realpart sig))))
	      (setf q (/ q (abs q)))
	      (format s "~6,3f ~6,3f~%"  i (let ((v (realpart sig)))
					     (cond ((< v -10) -100000)
						   ((<= v 10) v)
						   ((< 10 v) 100000)))
		     ))))
     #+inl(progn ;; rds
       (terpri s)
       (loop for pil across d and sig across f and i from 0 below 4000 do
	    (let ((q (complex (realpart pil)
			      (realpart sig))))
	      (setf q (/ q (abs q)))
	      (format s "~6,3f ~6,3f~%"  i (realpart sig)
		     )))))))


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

(defun realpartv (a)
  (let* ((n (length a))
	 (r (make-array n :element-type 'double-float)))
    (dotimes (i n)
      (setf (aref r i) (realpart (aref a i))))
    r))

(store-dfloat "/dev/shm/pilot3.dfloat" 
	      (realpartv (napa-fft:fft (fftshift *pilot3*))))
(store-dfloat "/dev/shm/rds.dfloat" 
	      (realpartv (napa-fft:fft (fftshift *rds*))))


#+nil
(let* ((n  (* 10 1024)) ;; verify data export with a sine wave
       (a (make-array n :element-type 'single-float)))
  (dotimes (i n)
    (setf (aref a i) (sin (* 2 #.(coerce pi 'single-float) i 200 (/ 1f0 n)))))
  (store-sfloat "/dev/shm/sin.sfloat" a)
  nil)

(time
 )

