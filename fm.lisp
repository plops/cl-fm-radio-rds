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

(defparameter *rate* 0)
(defparameter *n-complex* 0)
(defparameter *input* 0)

(progn
  (setf *rate* 2048000)
  (setf *n-complex* (floor (expt 2 22) #+nil 47448064 2))
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

(defparameter *spurious-freq*
  '((-130718.8 900) ;; second number is bandwidth (diameter)
    (-120718.8  500)
    (-119250.0  500)
    (-69625.0  530)
    (-59625.0  500)
    (-58187.5  500)
    (-16250.0  500)
    (52500.0  500)
    (62500.0  400)
    (63968.8  500)
    (82843.8  500)
    (113593.8  900)
    (123593.8  300)
    (125031.2  500)))


(defparameter *kin-filt*
  (let* ((n (length *kin*))
	 (nh (floor n 2))
	 
	 (bw (floor (* n 160000) ;; bandpass +/- 160kHz
		    *rate*))
	 (n-small (next-power-of-two (* 2 bw))
	   )
	 (bww (floor n-small 2))
	 (res (make-array n-small :element-type (array-element-type *kin*)
			  :initial-element (complex .0d0))))
    (format t "2*bw=~6,2f bww=~d n-small=~d~%" bw bww n-small)
    (loop for i from (+ (- bww) nh)
       below (+ bww nh) and ii from 0
       do
	 (setf (aref res ii) (if (< (- nh bw) i (+ nh bw))
				 (aref *kin* i)
				 (complex 0d0))))
    (loop for (c bw) in *spurious-freq* do
	 (let ((bw-bin (floor (* 2 n-small bw) ;; radius but i have to increase it
			      (* 2 *rate*)))
	       (c-bin (floor (* n c)
			     *rate*)))
	   
	   (let ((pos (+ (floor n-small 2)
			 c-bin)))
	     (loop for i from (- pos bw-bin)
	       upto (+ pos bw-bin) do
		  (when (<= 0 i (1- n-small))
		   (setf (aref res i) (complex 0d0)))))))
    res))

(progn ;; reverse ft
  (defparameter *input-filt*
    (napa-fft:ifft
     (fftshift *kin-filt*)))
  (store-cdfloat "/dev/shm/kin-filt.cdfload"
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
			  :element-type 'double-float)))
      (loop for i from 1 below n do
	   (let* ((s (aref in i))
		  (ds/dt (- s (aref in (1- i)))))
	     (declare (type (complex double-float) s ds/dt))
	     (setf (aref d i) (imagpart (/ ds/dt
					   s)))))
      d))
  (store-dfloat "/dev/shm/demod-heterodyn.dfloat" *demod-heterodyn*)
  nil)

(progn ;; cut out +/-50Hz around 19kHz
  (let* ((bw 100)
	 (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	 (nn (length d))
	 (center-bin (floor (* 19d3 nn
			       (/ 512d3))))
	 (band-bin (next-power-of-two
		    (* bw nn
		       (/ 512d3))))
	 (nh (floor band-bin 2))
	 (a (make-array (floor band-bin) :element-type (array-element-type d))))
    (loop for i from (- center-bin nh) below (+ center-bin nh) and ii from 0
       do
	 (setf (aref a ii) (aref d (+ (floor nn 2) i))))
    (store-cdfloat 
     (format nil "/dev/shm/pilot~6,3f.cdfloat" (* band-bin 512d3 (/ nn))) 
     (napa-fft:ifft (fftshift a)))
    nil))

(progn ;; cut out +/-2000Hz around 57kHz
  (let* ((bw 4000)
	 (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	 (nn (length d))
	 (center-bin (floor (* 57d3 nn
			       (/ 512d3))))
	 (band-bin (next-power-of-two
		    (* bw nn
		       (/ 512d3))))
	 (nh (floor band-bin 2))
	 (a (make-array (floor band-bin) :element-type (array-element-type d))))
    (loop for i from (- center-bin nh) below (+ center-bin nh) and ii from 0
       do
	 (setf (aref a ii) (aref d (+ (floor nn 2) i))))
    (store-cdfloat 
     (format nil "/dev/shm/rds~6,3f.cdfloat" (* band-bin 512d3 (/ nn))) 
     (napa-fft:ifft (fftshift a)))
    nil))

(defun cdf->df (a)
  (let ((b (make-array (length a) :element-type 'double-float)))
    (dotimes (i (length a))
      (setf (aref b i) (realpart (aref a i))))
    b))

(defparameter *pilot*
   (progn ;; cut out +/-50Hz around 19kHz, but leave at 19kHz
     (let* ((bw 100)
	    (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	    (nn (length d))
	    (nh (floor nn 2))
	    (center-bin (floor (* 19d3 nn
				  (/ 512d3))))
	    (band-bin (floor (* bw nn
				(/ 512d3))
			     2))
	    (a (make-array nn :element-type (array-element-type d))))
       (loop for i from (- center-bin band-bin)
	  below (+ center-bin band-bin)
	  do
	  (setf (aref a (+ nh i)) (aref d (+ nh i))
		(aref a (- nh i)) (aref d (- nh i))))
       (cdf->df (napa-fft:ifft (fftshift a))))))

(progn
 (defparameter *pilot-c* 
   (progn ;; cut out +/-50Hz around 19kHz, leave at 19kHz but single side band
     (let* ((bw 100)
	    (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	    (nn (length d))
	    (nh (floor nn 2))
	    (center-bin (floor (* 19d3 nn
				  (/ 512d3))))
	    (band-bin (floor (* bw nn
				(/ 512d3))
			     2))
	    (a (make-array nh :element-type (array-element-type d))))
       (loop for i from (- center-bin band-bin)
	  below (+ center-bin band-bin) and ii from 0
	  do
	  (setf (aref a (+ (floor nh 2) center-bin ii)) (aref d (+ nh i))))
       (napa-fft:ifft (fftshift a)))))
 (store-cdfloat "/dev/shm/pilot.cdfloat"
		*pilot-c*))

(defmacro with-plot ((stream fn) &body body)
  `(progn
     (with-open-file (s "/dev/shm/o.gp" :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
       (format s "plot ~s u 1:2 w l;pause -1" ,fn))
     (with-open-file (,stream ,fn :direction :output
		       :if-exists :supersede
		       :if-does-not-exist :create)
       ,@body)))

(let ((old-phi 0d0)
      (old-filt 0d0))
  (defun reset-dpll ()
    (setf old-phi 0d0
	  old-filt 0d0))
 (defun dpll (z)
   (declare (type (complex double-float) z)
	    (values (complex double-float) &optional))
   (let* ((f0 4800)
	  (fs (* 8 f0))
	  (eta .707)
	  (fn 10)
	  (omega_n (* 2 pi fn))
	  (c2 (* 2 eta omega_n (/ fs)))
	  (c1 (/ (expt c2 2)
		 (* 4 (expt eta 2))))
	  (c .1))
     (unless (and (< 0 c1)
		  (< (- (* 2 c2) 4) c1 c2))
       (error "filter parameters are not stable"))
     (let* ((phi_i (phase z)) ;; PD
	    (phi_e (- phi_i old-phi))) 
       (progn ;; digital filter
	 (let* ((top (+ (* c1 phi_e) old-filt))
		(bottom (* c2 phi_e))
		(filt-out (+ top bottom)))
	   (setf old-filt top)
	   (progn ;; VCO
	     (setf old-phi (+ c filt-out old-phi))
	     (exp (complex 0 old-phi)))))))))


#+nil
(with-plot (s "/dev/shm/o.dat")
  (reset-dpll)
  (loop for e across *pilot-c* and i below 30000 do
       (format s "~f ~9,4f~%" i (realpart (dpll e)))))


(defparameter *rds*
 (progn ;; cut out +/-2kHz around 57kHz, but leave at 57 kHz
   (let* ((bw 4000)
	  (d (fftshift (napa-fft:fft *demod-heterodyn*)))
	  (nn (length d))
	  (nh (floor nn 2))
	  (center-bin (floor (* 57d3 nn
				(/ 512d3))))
	  (band-bin (floor (* bw nn
			      (/ 512d3))
			     2))
	  (a (make-array nn :element-type (array-element-type d))))
     (loop for i from (- center-bin band-bin)
	below (+ center-bin band-bin)
	do
	  (setf (aref a (+ nh i)) (aref d (+ nh i))
		(aref a (- nh i)) (aref d (- nh i))))
       a)))



(progn ;; store both signals in file
  (store-dfloat "/dev/shm/rds.dfloat"
		(cdf->df (napa-fft:ifft (fftshift *rds*))))
  (store-dfloat "/dev/shm/pilot.dfloat"
		*pilot*))



#+nil
(sb-ext:gc :full t)

