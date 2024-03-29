(eval-when (:compile-toplevel :execute)
 (ql:quickload "napa-fft3"))

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
     (write-sequence o s)))
  (values))

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
     (write-sequence o s)))
  (values))

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
     (write-sequence o s)))
  (values))

(defun cdf->df (a)
  (let ((b (make-array (length a) :element-type 'double-float)))
    (dotimes (i (length a))
      (setf (aref b i) (realpart (aref a i))))
    b))

(defparameter *rate* 0)
(defparameter *n-complex* 0)
(defparameter *input* 0)

(progn ;; read in file
  (setf *rate* 2048000)
  (setf *n-complex* (floor (expt 2 26)))
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
	  (setf (aref c i) 
		(complex (- x 127d0)
			 (- y 127d0)))))
      c)))

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


(defparameter *kin*
  (fftshift
   (napa-fft:fft *input*)))

(defparameter *small-rate* 0)

(defparameter *kin-filt*
  (let* ((n (length *kin*))
	 (nh (floor n 2))
	 (center-bin (/ (* n -507250)
			*rate*))
	 (bw (floor (* n 180000) ;; bandpass +/- 90kHz
		    *rate*))
	 (small-n (* 2
		     1 ;; maybe oversample for later pll application
		     (next-power-of-two bw)))
	 (res (make-array small-n :element-type (array-element-type *kin*)
			  :initial-element (complex .0d0))))
    (defparameter *small-rate* (/ (* small-n *rate*)
				  n))
    (format t "smaller rate ~12,4f Hz~%" *small-rate*)
    (loop for i
       from  (- bw)
       below bw 
       and ii from (- (floor small-n 2)
		      bw)
       do
	 (setf (aref res ii) (aref *kin* (+ i nh center-bin))))
    res))

(progn ;; reverse ft
  (defparameter *input-filt*
    (napa-fft:ifft
     (fftshift *kin-filt*)))
  (store-cdfloat "/dev/shm/kin-filt.cdfloat"
		 *input-filt*)
  nil)



;; s = A e^ip
;; ds/dt = d/dt A e^ip = A ip' e^ip
;; ds/dt /s = ip'
;; p' = Im[ds/dt /s]
(progn ;; demodulate using heterodyne division
  (defparameter *demod-heterodyn*
    (let* ((in *input-filt*)
	   (n (length in))
	   (d (make-array n 
			  :element-type 'double-float))
	   (old-sample (complex 0d0)))
      (loop for i from 1 below n do
	   (let* ((s (aref in i))
		  (ds/dt (- s (aref in (1- i)))))
	     (declare (type (complex double-float) s ds/dt))
	     (setf old-sample
		   (setf (aref d i) (imagpart (if (< (abs s) .01d0)
						  old-sample
						  (/ ds/dt
						     s)))))))
      d))
  (store-dfloat "/dev/shm/demod-heterodyn.dfloat" *demod-heterodyn*)
  nil)

(progn
  (defparameter *pilot-c* 
    (progn ;; cut out +/-50Hz around 19kHz, leave at 19kHz but single side band
      (let* ((bw 100)
	     (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	     (nn (length d))
	     (nh (floor nn 2))
	     (center-bin (floor (* 19d3 nn
				   (/ *small-rate*))))
	     (band-bin (floor (* bw nn
				 (/ *small-rate*))
			      2))
	     (a (make-array nh :element-type '(complex double-float))))
	(loop for i from (- center-bin band-bin)
	   below (+ center-bin band-bin) and ii from (- band-bin)
	   do
	     (setf (aref a (+ (floor nh 2) center-bin ii)) (aref d (+ nh i))))
	(napa-fft:ifft (fftshift a)))))
  (store-cdfloat "/dev/shm/pilot.cdfloat"
		 *pilot-c*))

(progn
  (defparameter *rds-c* 
    (progn ;; cut out +/-2000Hz around 57kHz, leave at 57kHz but single side band
      (let* ((bw 4000)
	     (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	     (nn (length d))
	     (nh (floor nn 2))
	     (center-bin (floor (* 57d3 nn
				   (/ *small-rate*))))
	     (band-bin (floor (* bw nn
				 (/ *small-rate*))
			      2))
	     (a (make-array nh :element-type '(complex double-float))))
	(loop for i from (- center-bin band-bin)
	   below (+ center-bin band-bin) and ii from (- band-bin)
	   do
	     (setf (aref a (+ (floor nh 2) center-bin ii)) (aref d (+ nh i))))
	(napa-fft:ifft (fftshift a)))))
  (store-cdfloat "/dev/shm/rds.cdfloat"
		 *pilot-c*))

(sb-ext:gc :full t)

(let ((old-phi 0d0)
      (old-phi3-cont 0d0)
      (old-phi3 0d0)
      (old-filt 0d0)
      (divider/ 1d0)
      (c2 0d0)
      (c1 0d0)
      (c 0d0))
  (declare (type double-float old-phi old-phi3 
		 old-phi3-cont old-filt c2 c1 
		 c divider/))
  (defun reset-dpll3 (&key (z (complex 0d0)) (z3 (complex 0d0))
		      (f0 -19.0d3)
		      (divider 1d0)
		      (fs 256d3) (eta .707d0) (fn 5000d0))
    (setf old-phi (phase z)
	  old-phi3 (phase z3)
	  old-phi3-cont old-phi3
	  old-filt 0d0
	  divider/ (/ divider)
	  c2 (* 2 eta (* 2 pi fn) (/ fs))
	  c1 (/ (expt c2 2)
		(* 4 (expt eta 2)))
	  c (/ (* 2 pi (* divider f0)) 
	       fs))
    (unless (and (< 0 c1)
		 (< (- (* 2 c2) 4) c1 c2))
      (error "filter parameters are not stable")))
  (defun dpll3 (z)
    (declare (type (complex double-float) z)
	     (optimize speed)
	     (values double-float 
		     double-float
		     (complex double-float) 
		     (complex double-float) &optional))
    (let* ((phi_i (phase z)) ;; PD
	   (phi_e (- phi_i old-phi)))
      
      (if (< phi_e (* .9 -2d0 pi)) ;; phase unwrapping
	  (incf phi_e (* 2 pi))
	  (if (<  (* .9 2d0 pi) phi_e)
	      (decf phi_e (* 2 pi))))
      
					; (format t "~a~%" phi_e)
      (progn				;; digital filter
	(let* ((top (+ (* c1 phi_e) old-filt))
	       (bottom (* c2 phi_e))
	       (filt-out (+ top bottom)))
	  (setf old-filt top)
	  (let* ((znew	;; VCO
		  (exp (complex 0 (+ (* divider/ c) filt-out old-phi))))
		 (znew3 (exp (complex 0
				      (+ c filt-out old-phi3)))))
	    (setf old-phi (phase znew)
		  old-phi3-cont (+ c filt-out old-phi3-cont)
		  old-phi3 (phase znew3))
	    (values
	     old-phi3-cont
	     filt-out
	     znew
	     znew3)))))))

(progn ;; calculate phase derivative
  (let* ((n (length *pilot-c*))
	 (ac (make-array n :element-type '(complex double-float)))
	 (a (make-array n :element-type 'double-float))
	 (ph (make-array n :element-type 'double-float))
	 (max-err 0d0))
    (reset-dpll3 :z (aref *pilot-c* 0)
		 :z3 (aref *rds-c* 0)
		 :f0 (* 2 -19d3)
		 :divider 3d0
		 :fs *small-rate*
		 :fn 20d0)
    (dotimes (i (1- n))
      (let ((e (aref *pilot-c* i))
	    (r-next (aref *rds-c* (1+ i))))
	(multiple-value-bind (phi er z z3) (dpll3 e)
	  (setf (aref ac i) (/ r-next z3)
		(aref ph i) phi
		max-err (max max-err (abs er))))))
    (dotimes (i (1- n))
      (setf (aref a i) (imagpart (/ (- (aref ac (1+ i))
				       (aref ac i))
				    (aref ac i)))))
    (defparameter *dphase* a)
    (defparameter *dphase-arg* ph)
    (defparameter *dphase-max-err* max-err)
    (store-dfloat "/dev/shm/dphase.dfloat" a)))

(defmacro with-plot ((stream fn) &body body)
  `(progn
     (with-open-file (s "/dev/shm/o.gp" :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
       (format s #+nil "plot ~s u 1:2 w l, ~s u 1:3 w l, ~s u 1:4 w l, ~s u 1:5 w l;pause -1"
	       #+nil "plot ~s u 1:2 w p pt 0 ps 0;pause -1"
	       "plot ~s u 1:2 w l;pause -1"
	       ,fn))
     (with-open-file (,stream ,fn :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
       ,@body)))

#+nil
(with-plot (s "/dev/shm/o.dat")
  (loop for i from 10000 below 30000 do
       (let ((e (aref *dphase* i))
	     (x (aref *dphase-arg* i)))
	 (format s "~f ~f~%" 
		 (let ((s (* 2 pi 48))) 
		   (- (/ x s)) ) 
		 e))))




