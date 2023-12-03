import openseespy.opensees as ops
def applyGM(startTime):
	pass
	# GMfile='./fromUser/GroundMotion.acc'
	# # Uniform EXCITATION: acceleration input
	# IDloadTag = 400			# load tag
	# # dt = 0.00001			# time step for input ground motion
	# # maxNumIter = 10
	# # GMdirection=1
	# # Tol=1e-3

	# TSTAG=2
	# period=0.01
	# factor=5e-4
	# tStart=0.01
	# tEnd=0.25
	# #ops.timeSeries('Path', TSTAG, '-dt', dt, '-filePath', GMfile, '-factor', GMfact, '-useLast', '-startTime', startTime)
	# ops.timeSeries('Triangle', TSTAG, tStart, tEnd, period, '-factor', factor)
	# #ops.pattern('UniformExcitation', IDloadTag, GMdirection, '-accel', 2) 

	# ops.pattern('MultipleSupport', IDloadTag-1)

	# ops.groundMotion(IDloadTag, 'Plain', '-disp',TSTAG, '-int', 'Trapezoidal', '-fact', 1.0) 
	# ops.imposedMotion(0, 1,IDloadTag)
	
def removeGM():
	pass
	# TSTAG=2
	# IDloadTag = 400
	# ops.remove('loadPattern',IDloadTag-1) 
	# ops.remove('loadPattern',IDloadTag)
	# ops.remove('timeSeries', TSTAG) 