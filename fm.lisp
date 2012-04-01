
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
     (dotimes (i (floor n 2))
       (setf (aref d i) (atan (imagpart (aref c i))
			      (realpart (aref c i)))))
     d)))

(time
 (with-open-file (s "/dev/shm/demod.sfload"
		    :direction :output
		    :element-type '(unsigned-byte 8)
		    :if-does-not-exist :create
		    :if-exists :supersede)
   (let* ((a (sb-alien:sap-alien (sb-sys:vector-sap *demod*)
				 (sb-alien:* sb-alien:unsigned-char)))
	  (n (length *demod*))
	  (o (make-array (* n 4) :element-type '(unsigned-byte 8))))
     (dotimes (i (* 4 n))
       (setf (aref o i) (deref a i)))
     (write-sequence o s))))