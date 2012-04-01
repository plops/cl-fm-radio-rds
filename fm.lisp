
(time
 (defparameter *demod*
   (let* ((n 63438848)
	  (a (make-array n :element-type '(unsigned-byte 8)))
	  (c (make-array (floor n 2) :element-type '(complex single-float)))
	  (d (make-array (length c)
			 :element-type 'single-float)))
     (with-open-file (s "/home/martin/Downloads2/dlf.fm"
			:element-type '(unsigned-byte 8))
       (read-sequence a s))
     (dotimes (i (floor n 2))
       (setf (aref c i) 
	     (complex (- (aref a (* 2 i)) 127f0)
		      (- (aref a (1+ (* 2 i))) 127f0))))
     (loop for i from 1 below (floor n 2) do
	  (let ((prod (* (aref c i) (conjugate (aref c (1- i))))))
	   (setf (aref d i) (realpart prod) #+nil(atan (imagpart prod)
				  (realpart prod)))))
     d)))

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