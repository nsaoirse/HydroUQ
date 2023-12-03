# # # -*- coding: utf-8 -*-
# # """
# # Created on Wed Sep 13 13:54:20 2023

# # @author: nsaoirse
# # """


# import openseespy.opensees as ops
# from openseespy.opensees import *
# import numpy as np 

# class FOAMySeesInstance():
	# def __init__(self, iNPUTtIMESTEP, config, parent=None):
		# pass
ops.wipe()


beamType='elastic'



beamType='dispBeamCol'


beamType='forceBeamCol'
beamType='NLBeamCol'


noBraces=0
useBeams=1

constraintType='RIGIDLINK'	
constraintType='EQDOFMIXED'
constraintType='EQDOF'



useSlabs=1

BaseIsolated=0

zBumpDebug=0.0

numDOF=6
rcdofs=[1,1,2,2,3,3,4,4,5,5,6,6]
dofs=[1,2,3,4,5,6]
FOAMySeesInstance.osi=ops.model('basic','-ndm',3,'-ndf',numDOF)


FOAMySeesInstance.coupledNodes=[]

nElemBOTTOM=2
nElem2nd=2
nElemTOP=2

numX=4
numY=4

## Materials
matTag=600
Fy=56*6895000 #ksi*conversion
E0=29000*6895000 #ksi*conversion
b=0.01
params=[18.00000, 0.92500, 0.15000]
ops.uniaxialMaterial('Steel02', matTag, Fy, E0, b, *params)
matTag=601
Fy=50*6895000 #ksi*conversion
E0=29000*6895000 #ksi*conversion
b=0.01
params=[18.00000, 0.92500, 0.15000]
ops.uniaxialMaterial('Steel02', matTag, Fy, E0, b, *params)

matTag=650
fpc=-7.2*6895000 #ksi*conversion
epsc0=-0.00326
fpcu=-1.44000*6895000 #ksi*conversion
epsU=-0.01631
lambda1=0.10000
ft=0.63640*6895000 #ksi*conversion
Ets=290.47375*6895000 #ksi*conversion

ops.uniaxialMaterial('Concrete02', matTag, fpc, epsc0, fpcu, epsU, lambda1, ft, Ets)
OuterDiam=4
CFTWallT=0.5

## Sections]
secTag=501
ops.section('Fiber', secTag, '-GJ', 10000000000.00000)

matTag=650
numSubdivCirc=8
numSubdivRad=8
center=[0.00000, 0.00000]
rad=[0.00000, (OuterDiam-CFTWallT)*0.0254] #INCH*conversionToMeters
ang=[0.00000, 360.00000]

ops.patch('circ', matTag, numSubdivCirc, numSubdivRad, *center, *rad, *ang)

matTag=600
numSubdivCirc=8
numSubdivRad=4
center=[0.00000, 0.00000]
rad=[(OuterDiam-CFTWallT)*0.0254, OuterDiam*0.0254] #INCH*conversionToMeters
ang=[0.00000, 360.00000]
ops.patch('circ', matTag, numSubdivCirc, numSubdivRad, *center, *rad, *ang)

############################

ops.beamIntegration('Lobatto', 1, secTag, 2)

# FIRST STORY COLUMNS
##################################################################################################################################################################
pileMass=7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519

slabMass=1000

nElem=nElemBOTTOM
Zbase=2.0
Ztop=2.6180000
Xmin=40.88230
Xmax=41.89830
Ymin=-0.5080000
Ymax=0.5080000


BeamConnLength=(Ymax-Ymin)/numX

pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]



ConnectionLength=0.05
pc=0
for location in pileLocs:

	node1=[location[0], location[1],  Zbase]
	node2=[location[0], location[1], Ztop-ConnectionLength]

	beamNormal=[-1,0,0]


	xNodeList=np.linspace(node1[0],node2[0],nElem+1)
	yNodeList=np.linspace(node1[1],node2[1],nElem+1)
	zNodeList=np.linspace(node1[2],node2[2],nElem+1)
	nodalMass=(pileMass)/len(xNodeList)

	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
		FOAMySeesInstance.coupledNodes.append(nodeNum)

	coordTransf = 'Corotational'

	coordTransf='Linear'
	coordTransf='PDelta'
	#############################

	## Model



	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, A, E, G, J, Iy, Iz, nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1

if useSlabs==1:
	# FIRST STORY SLAB
	##################################################################################################################################################################

	zBumpDebug=0.0
	Ztop+=zBumpDebug


	pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]


	eleType='shell'
	eleArgs=[100]

	matTag=625

	h=0.0254*1/2	 #SLAB THICKNESS IN METERS> WHAT IS THIS? CHECK KENS THESIS

	nu=0.3

	E=8e9 #plywood panel

	E=2e11 #steel panel

	rho=1000   #SLAB DENSITY > WHAT IS THIS? CHECK KENS THESIS

	ops.nDMaterial('ElasticIsotropic', matTag, E, nu)

	secTag=100

	ops.section('ElasticMembranePlateSection', secTag, E, nu, h, rho)
	# ops.node(51, pileLocs[0][0],pileLocs[0][1],pileLocs[0][2])
	# ops.node(1051,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2])
	# ops.node(2051,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2])
	# ops.node(3051,pileLocs[3][0],pileLocs[3][1],pileLocs[3][2])



	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]



	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=100000
	startEle=100000

	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)
	SlabNodes=np.linspace(startNode,startNode+((numX+1)*(numY+1))-1,((numX+1)*(numY+1)))
	for nodeNum in SlabNodes:
		FOAMySeesInstance.coupledNodes.append(int(nodeNum))


Ztop-=zBumpDebug

integTag=15
	
ops.beamIntegration('Lobatto', integTag, secTag, 2)

if useBeams==1:
	
	# 1st STORY BEAMS
	##################################################################################################################################################################
	
	
	pileLocs=[[[Xmin,Ymin+BeamConnLength,Ztop],[Xmin,Ymax-BeamConnLength,Ztop]],[[Xmin+BeamConnLength,Ymax,Ztop],[Xmax-BeamConnLength,Ymax,Ztop]],[[Xmax,Ymax-BeamConnLength,Ztop],[Xmax,Ymin+BeamConnLength,Ztop]],[[Xmax-BeamConnLength,Ymin,Ztop],[Xmin+BeamConnLength,Ymin,Ztop]]]
	
	
	
	pileMass=7840*0.84*0.00064516*(Xmax-Xmin)
	pc=10
	
	# Nominal Size 3)	Weight	Wall Thickness   Area	 I 	J
	# (in x in x in)	(lbf/ft)	(in)		(in2)	(in4)  (in4)
	#  2 x 2 x 1/8		 3.05		0.116	  	0.84*0.00064516	0.486*0.00064516*0.00064516	0.796*0.00064516*0.00064516
	
	# section('Elastic', secTag, E_mod, A, Iz, G_mod=None)
	A=0.84*0.00064516
	Iz=0.486*0.00064516*0.00064516
	Iy=0.486*0.00064516*0.00064516
	Jxx=0.796*0.00064516*0.00064516
	E_mod=29000*6895000 #ksi*conversion
	G_mod=E_mod/(2*(1.3))
	secTag=15
	ops.section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)

	for location in pileLocs:
	
		node1=[location[0][0], location[0][1],  location[0][2]]
		node2=[location[1][0], location[1][1],  location[1][2]]
	
		beamNormal=[0,0,1]
	
	
	
	
		xNodeList=np.linspace(node1[0],node2[0],numX-1)
		yNodeList=np.linspace(node1[1],node2[1],numX-1)
		zNodeList=np.linspace(node1[2],node2[2],numX-1)
		nodalMass=(pileMass)/len(xNodeList)
	
		for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
			ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
	
		coordTransf = 'Corotational'
	
		coordTransf='Linear'
		coordTransf='PDelta'
		#############################
	
		## Model
	
	
	
		nodRotMass=0.
		for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
			ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])
	
		for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
			ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
			if beamType=='elastic':
				ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, A, E, G, J, Iy, Iz, nodeNum+100000)
			elif beamType=='dispBeamCol':
				ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
			elif beamType=='forceBeamCol':
				ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
			elif beamType=='NLBeamCol':
				ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)	   
		pc+=1
	




# 2nd STORY COLUMNS
##################################################################################################################################################################
pileMass=7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(18*0.0254)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(18*0.0254)*3.141519




Zbase=2.6180000
#Ztop=1.8288
Ztop=2.6180000+(18*0.0254)
Xmin=40.88230
Xmax=41.89830
Ymin=-0.5080000
Ymax=0.5080000

pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]
ConnectionLength=0.05
pc=5
for location in pileLocs:

	node1=[location[0], location[1],  Zbase]
	node2=[location[0], location[1], Ztop-ConnectionLength]

	beamNormal=[-1,0,0]


	xNodeList=np.linspace(node1[0],node2[0],nElem2nd+1)
	yNodeList=np.linspace(node1[1],node2[1],nElem2nd+1)
	zNodeList=np.linspace(node1[2],node2[2],nElem2nd+1)
	nodalMass=(pileMass)/len(xNodeList)

	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
		FOAMySeesInstance.coupledNodes.append(nodeNum)

	coordTransf = 'Corotational'

	coordTransf='Linear'
	coordTransf='PDelta'
	#############################

	## Model



	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, A, E, G, J, Iy, Iz, nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1







if useSlabs==1:
	# 2nd STORY SLAB
	##################################################################################################################################################################

	zBumpDebug=0.0
	Ztop+=zBumpDebug
	pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]
	eleType='shell'
	eleArgs=[100]

	# ops.node(51, pileLocs[0][0],pileLocs[0][1],pileLocs[0][2])
	# ops.node(1051,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2])
	# ops.node(2051,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2])
	# ops.node(3051,pileLocs[3][0],pileLocs[3][1],pileLocs[3][2])



	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]


	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=200000
	startEle=200000

	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)

Ztop-=zBumpDebug


if useBeams==1:
	
	# 2nd STORY BEAMS
	##################################################################################################################################################################
	
	
	pileLocs=[[[Xmin,Ymin+BeamConnLength,Ztop],[Xmin,Ymax-BeamConnLength,Ztop]],[[Xmin+BeamConnLength,Ymax,Ztop],[Xmax-BeamConnLength,Ymax,Ztop]],[[Xmax,Ymax-BeamConnLength,Ztop],[Xmax,Ymin+BeamConnLength,Ztop]],[[Xmax-BeamConnLength,Ymin,Ztop],[Xmin+BeamConnLength,Ymin,Ztop]]]
	
	pileMass=7840*0.84*0.00064516*(Xmax-Xmin)
	pc=20
	
	# Nominal Size 3)	Weight	Wall Thickness   Area	 I 	J
	# (in x in x in)	(lbf/ft)	(in)		(in2)	(in4)  (in4)
	#  2 x 2 x 1/8		 3.05		0.116	  	0.84*0.00064516	0.486*0.00064516*0.00064516	0.796*0.00064516*0.00064516
	
	# section('Elastic', secTag, E_mod, A, Iz, G_mod=None)
	A=0.84*0.00064516
	Iz=0.486*0.00064516*0.00064516
	Iy=0.486*0.00064516*0.00064516
	Jxx=0.796*0.00064516*0.00064516
	E_mod=29000*6895000 #ksi*conversion
	G_mod=E_mod/(2*(1.3))
	secTag=15
	#ops.section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)
	
	#ops.beamIntegration('Lobatto', 15, secTag, 2)
	for location in pileLocs:
	
		node1=[location[0][0], location[0][1],  location[0][2]]
		node2=[location[1][0], location[1][1],  location[1][2]]
	
		beamNormal=[0,0,1]
	
	
		xNodeList=np.linspace(node1[0],node2[0],numX-1)
		yNodeList=np.linspace(node1[1],node2[1],numX-1)
		zNodeList=np.linspace(node1[2],node2[2],numX-1)
		nodalMass=(pileMass)/len(xNodeList)
	
		for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
			ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
	
		coordTransf = 'Corotational'
	
		coordTransf='Linear'
		coordTransf='PDelta'
		#############################
	
		## Model
	
	
	
		nodRotMass=0.
		for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
			ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])
	
		for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
			ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
			if beamType=='elastic':
				ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, A, E, G, J, Iy, Iz, nodeNum+100000)
			elif beamType=='dispBeamCol':
				ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
			elif beamType=='forceBeamCol':
				ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
			elif beamType=='NLBeamCol':
				ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag, pc*1000+nodeNum+100000,'-iter', 1)
		pc+=1
	



# TOP STORY COLUMNS
##################################################################################################################################################################
pileMass=7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(18*0.0254)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(18*0.0254)*3.141519

Zbase=2.6180000+(18*0.0254)
#Ztop=1.8288
Ztop=2.6180000+(18*0.0254)*2
Xmin=40.88230
Xmax=41.89830
Ymin=-0.5080000
Ymax=0.5080000



pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]

pc=15
for location in pileLocs:

	node1=[location[0], location[1],  Zbase]
	node2=[location[0], location[1], Ztop]

	beamNormal=[-1,0,0]


	xNodeList=np.linspace(node1[0],node2[0],nElemTOP+1)
	yNodeList=np.linspace(node1[1],node2[1],nElemTOP+1)
	zNodeList=np.linspace(node1[2],node2[2],nElemTOP+1)
	nodalMass=(pileMass)/len(xNodeList)

	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
		FOAMySeesInstance.coupledNodes.append(nodeNum)

	coordTransf = 'Corotational'

	coordTransf='Linear'
	coordTransf='PDelta'
	#############################

	## Model



	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, A, E, G, J, Iy, Iz, nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag, pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1


if useSlabs==1:
	# TOP STORY SLAB
	##################################################################################################################################################################

	zBumpDebug=0.0
	Ztop+=zBumpDebug

	pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]


	eleType='shell'
	eleArgs=[100]

	# ops.node(51, pileLocs[0][0],pileLocs[0][1],pileLocs[0][2])
	# ops.node(1051,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2])
	# ops.node(2051,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2])
	# ops.node(3051,pileLocs[3][0],pileLocs[3][1],pileLocs[3][2])



	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]



	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=300000
	startEle=300000

	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)


Ztop-=zBumpDebug

if useBeams==1:

	# TOP STORY BEAMS
	##################################################################################################################################################################
	
	
	
	pileLocs=[[[Xmin,Ymin+BeamConnLength,Ztop],[Xmin,Ymax-BeamConnLength,Ztop]],[[Xmin+BeamConnLength,Ymax,Ztop],[Xmax-BeamConnLength,Ymax,Ztop]],[[Xmax,Ymax-BeamConnLength,Ztop],[Xmax,Ymin+BeamConnLength,Ztop]],[[Xmax-BeamConnLength,Ymin,Ztop],[Xmin+BeamConnLength,Ymin,Ztop]]]
	
	pileMass=7840*0.84*0.00064516*(Xmax-Xmin)
	
	
	pc=30
	
	# Nominal Size 3)	Weight	Wall Thickness   Area	 I 	J
	# (in x in x in)	(lbf/ft)	(in)		(in2)	(in4)  (in4)
	#  2 x 2 x 1/8		 3.05		0.116	  	0.84*0.00064516	0.486*0.00064516*0.00064516	0.796*0.00064516*0.00064516
	
	# section('Elastic', secTag, E_mod, A, Iz, G_mod=None)
	A=0.84*0.00064516
	Iz=0.486*0.00064516*0.00064516
	Iy=0.486*0.00064516*0.00064516
	Jxx=0.796*0.00064516*0.00064516
	E_mod=29000*6895000 #ksi*conversion
	G_mod=E_mod/(2*(1.3))
	secTag=15
	#ops.section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)
	
	#ops.beamIntegration('Lobatto', 15, secTag, 2)
	for location in pileLocs:
	
		node1=[location[0][0], location[0][1],  location[0][2]]
		node2=[location[1][0], location[1][1],  location[1][2]]
	
		beamNormal=[0,0,1]
	
		BeamConnLength=max(abs(location[0][0]-location[1][0]),abs(location[0][1]-location[1][1]))/numX
	
		xNodeList=np.linspace(node1[0],node2[0],numX-1)
		yNodeList=np.linspace(node1[1],node2[1],numX-1)
		zNodeList=np.linspace(node1[2],node2[2],numX-1)
		nodalMass=(pileMass)/len(xNodeList)
	
		for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
			ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
	
		coordTransf = 'Corotational'
	
		coordTransf='Linear'
		coordTransf='PDelta'
		#############################
	
		## Model
	
	
	
		nodRotMass=0.
		for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
			ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])
	
		for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
			ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
			if beamType=='elastic':
				ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, A, E, G, J, Iy, Iz, nodeNum+100000)
			elif beamType=='dispBeamCol':
				ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
			elif beamType=='forceBeamCol':
				ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
			elif beamType=='NLBeamCol':
				ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
		pc+=1


########### REPLACE WITH CONNECTION ELEMENTS############

#element('beamColumnJoint', eleTag, *eleNodes, Mat1Tag, Mat2Tag, Mat3Tag, Mat4Tag, Mat5Tag, Mat6Tag, Mat7Tag, Mat8Tag, Mat9Tag, Mat10Tag, Mat11Tag, Mat12Tag, Mat13Tag)
	# eleTag (int)	unique element object tag
	# eleNodes (list (int))	a list of four element nodes
	# Mat1Tag (int)	uniaxial material tag for left bar-slip spring at node 1
	# Mat2Tag (int)	uniaxial material tag for right bar-slip spring at node 1
	# Mat3Tag (int)	uniaxial material tag for interface-shear spring at node 1
	# Mat4Tag (int)	uniaxial material tag for lower bar-slip spring at node 2
	# Mat5Tag (int)	uniaxial material tag for upper bar-slip spring at node 2
	# Mat6Tag (int)	uniaxial material tag for interface-shear spring at node 2
	# Mat7Tag (int)	uniaxial material tag for left bar-slip spring at node 3
	# Mat8Tag (int)	uniaxial material tag for right bar-slip spring at node 3
	# Mat9Tag (int)	uniaxial material tag for interface-shear spring at node 3
	# Mat10Tag (int)	uniaxial material tag for lower bar-slip spring at node 4
	# Mat11Tag (int)	uniaxial material tag for upper bar-slip spring at node 4
	# Mat12Tag (int)	uniaxial material tag for interface-shear spring at node 4
	# Mat13Tag (int)	uniaxial material tag for shear-panel




numX-=2
# 1st Story Beams to Columns
BCNodeConns=[]

numX-=2
# 1st Story Beams to Columns
BCNodeConns=[[5000,10000],[5000,13000+numX],[6000,11000],[6000,10000+numX],[7000,12000],[7000,11000+numX],[8000,13000],[8000,12000+numX]]

# 2nd Story Beams to Columns
for x in [[15000,20000],[15000,23000+numX],[16000,21000],[16000,20000+numX],[17000,22000],[17000,21000+numX],[18000,23000],[18000,22000+numX]]:
	BCNodeConns.append(x)

# Top Story Beams to Columns
for x in [[15000+nElemTOP,30000],[15000+nElemTOP,33000+numX],[16000+nElemTOP,31000],[16000+nElemTOP,30000+numX],[17000+nElemTOP,32000],[17000+nElemTOP,31000+numX],[18000+nElemTOP,33000],[18000+nElemTOP,32000+numX]]:
	BCNodeConns.append(x)

# 1st to 2nd story Column to column connection

CCNodeConns=[[nElemBOTTOM,5000],[1000+nElemBOTTOM,6000],[2000+nElemBOTTOM,7000],[3000+nElemBOTTOM,8000]]


# 2nd to 3rd story Column to column connection

for x in [[nElem2nd+5000,15000],[nElem2nd+6000,16000],[nElem2nd+7000,17000],[nElem2nd+8000,18000]]:
	CCNodeConns.append(x)

direcs=[1,2,3,4,5,6]

matTags=[600,600,600,600,600,600]

eleTag=500000


if useBeams==1:
	
	for BCpair in BCNodeConns:
		rNodeTag=BCpair[0]
		cNodeTag=BCpair[1]
        
		#ops.element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *direcs)
		#ops.element('twoNodeLink', int(eleTag), *[BCpair[0],BCpair[1]], '-mat', *matTags, '-dir', *dirs,'-orient', *vecx, *vecyp)
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		
	
		#ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		eleTag+=1


for CCpair in CCNodeConns:
	rNodeTag=CCpair[0]
	cNodeTag=CCpair[1]
	eleNodes=[CCpair[0],CCpair[1]]

	ops.geomTransf(coordTransf, eleTag,-1,0,0)
	ops.element('forceBeamColumn', eleTag, *eleNodes, eleTag, 1)

	# ops.element('zeroLength', eleTag, *eleNodes, '-mat', *matTags, '-dir', *direcs)
	eleTag+=1


if useBeams==1:
	
	numX+=2
	dofs=[1,2,3,4,5,6]
	############### FIRST FLOOR SLAB TO BEAM CONNECTIONS
	AllSlabNodes=[]
	MNC=numX-1
	print(numX,numY)
	SN=np.array(np.linspace(100000,100000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))
	
	print(SN)
	
	SN=np.reshape(SN,((numX+1),(numY+1)))
	print(SN)
	
	LS1=list(SN[0,:])
	print(LS1)
	LS2=list(SN[:,numX])
	print(LS2)
	LS3=list(SN[numX,:])[::-1]
	print(LS3)
	LS4=list(SN[:,0])[::-1]
	print(LS4)
	
	AllSlabNodes.append(list(LS1[1:numX]))
	AllSlabNodes.append(list(LS2[1:numX]))
	AllSlabNodes.append(list(LS3[1:numX]))
	AllSlabNodes.append(list(LS4[1:numX]))
	
	with open('DEBUGMESH.txt', 'a+') as f:
	
		print(AllSlabNodes,file=f)
	
	LT1=np.linspace(10000,10000+numX-2,MNC)
	LT2=np.linspace(11000,11000+numX-2,MNC)
	LT3=np.linspace(12000,12000+numX-2,MNC)
	LT4=np.linspace(13000,13000+numX-2,MNC)
	
	AllTrussNodes=[]
	
	AllTrussNodes.append(list(LT1[:]))
	AllTrussNodes.append(list(LT2[:]))
	AllTrussNodes.append(list(LT3[:]))
	AllTrussNodes.append(list(LT4[:]))
	
	with open('DEBUGMESH.txt', 'a+') as f:
	
		print(AllTrussNodes,file=f)
	
	eqDOFS=[]
	for nodeGroup in range(0,4):
		for ind in range(0,MNC):
			coup=[int(AllSlabNodes[nodeGroup][ind]) ,int(AllTrussNodes[nodeGroup][ind])]
			eqDOFS.append(coup)
	
	dofs=[1,1,1,1,1,1]
	CT=6000000
	coordTransf='PDelta'
	for eqDOFpair in eqDOFS:
		rNodeTag=eqDOFpair[0]
		cNodeTag=eqDOFpair[1]
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		
		# beamNormal=[1,0,0]
		# ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
		# ops.element('forceBeamColumn', CT, *[rNodeTag, cNodeTag], CT, 15)
		#ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		# CT+=1


	
	############### SECOND FLOOR SLAB TO BEAM CONNECTIONS
	AllSlabNodes=[]
	
	print(numX,numY)
	SN=np.array(np.linspace(200000,200000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))
	
	print(SN)
	
	SN=np.reshape(SN,((numX+1),(numY+1)))
	print(SN)
	
	LS1=list(SN[0,:])
	print(LS1)
	LS2=list(SN[:,numX])
	print(LS2)
	LS3=list(SN[numX,:])[::-1]
	print(LS3)
	LS4=list(SN[:,0])[::-1]
	print(LS4)
	
	AllSlabNodes.append(list(LS1[1:numX]))
	AllSlabNodes.append(list(LS2[1:numX]))
	AllSlabNodes.append(list(LS3[1:numX]))
	AllSlabNodes.append(list(LS4[1:numX]))
	
	with open('DEBUGMESH.txt', 'a+') as f:
	
		print(AllSlabNodes,file=f)
	
	LT1=np.linspace(20000,20000+numX-2,MNC)
	LT2=np.linspace(21000,21000+numX-2,MNC)
	LT3=np.linspace(22000,22000+numX-2,MNC)
	LT4=np.linspace(23000,23000+numX-2,MNC)
	
	AllTrussNodes=[]
	
	AllTrussNodes.append(list(LT1[:]))
	AllTrussNodes.append(list(LT2[:]))
	AllTrussNodes.append(list(LT3[:]))
	AllTrussNodes.append(list(LT4[:]))
	
	with open('DEBUGMESH.txt', 'a+') as f:
	
		print(AllTrussNodes,file=f)
	eqDOFS=[]
	for nodeGroup in range(0,4):
		for ind in range(0,MNC):
			coup=[int(AllSlabNodes[nodeGroup][ind]) ,int(AllTrussNodes[nodeGroup][ind])]
			eqDOFS.append(coup)
	
	CT=5000000
	coordTransf='PDelta'
	for eqDOFpair in eqDOFS:
		rNodeTag=eqDOFpair[0]
		cNodeTag=eqDOFpair[1]
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		
		# beamNormal=[1,0,0]
		# ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
		# ops.element('forceBeamColumn', CT, *[rNodeTag, cNodeTag], CT, 15)
		#ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		# CT+=1
	
	
	
	############### THIRD FLOOR SLAB TO BEAM CONNECTIONS
	AllSlabNodes=[]
	
	print(numX,numY)
	SN=np.array(np.linspace(300000,300000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))
	
	print(SN)
	
	SN=np.reshape(SN,((numX+1),(numY+1)))
	print(SN)
	
	LS1=list(SN[0,:])
	print(LS1)
	LS2=list(SN[:,numX])
	print(LS2)
	LS3=list(SN[numX,:])[::-1]
	print(LS3)
	LS4=list(SN[:,0])[::-1]
	print(LS4)
	
	AllSlabNodes.append(list(LS1[1:numX]))
	AllSlabNodes.append(list(LS2[1:numX]))
	AllSlabNodes.append(list(LS3[1:numX]))
	AllSlabNodes.append(list(LS4[1:numX]))
	with open('DEBUGMESH.txt', 'a+') as f:
	
		print(AllSlabNodes,file=f)
	
	
	LT1=np.linspace(30000,30000+numX-2,MNC)
	LT2=np.linspace(31000,31000+numX-2,MNC)
	LT3=np.linspace(32000,32000+numX-2,MNC)
	LT4=np.linspace(33000,33000+numX-2,MNC)
	
	
	
	AllTrussNodes=[]
	
	AllTrussNodes.append(list(LT1[:]))
	AllTrussNodes.append(list(LT2[:]))
	AllTrussNodes.append(list(LT3[:]))
	AllTrussNodes.append(list(LT4[:]))
	
	with open('DEBUGMESH.txt', 'a+') as f:
	
		print(AllTrussNodes,file=f)
	eqDOFS=[]
	for nodeGroup in range(0,4):
		for ind in range(0,MNC):
			coup=[int(AllSlabNodes[nodeGroup][ind]) ,int(AllTrussNodes[nodeGroup][ind])]
			eqDOFS.append(coup)
	
	CT=4000000
	coordTransf='PDelta'
	for eqDOFpair in eqDOFS:
		rNodeTag=eqDOFpair[0]
		cNodeTag=eqDOFpair[1]
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		
		# beamNormal=[1,0,0]
		# ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
		# ops.element('forceBeamColumn', CT, *[rNodeTag, cNodeTag], CT, 15)
		#ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		# CT+=1
	

if noBraces==1:
	pass
else:
	############### Chevron braces #####################


	# corner nodes nElemBOTTOM  x 1000+nElemBOTTOM  y  2000+nElemBOTTOM  x  3000+nElemBOTTOM   y   nElemBOTTOM

	# mid nodes			   20000+numX/2	 21000+numX/2	 22000+numX/2		  23000+numX/2
			


	CT=101010101
	eleNo=CT
	beamNormal=[1,0,0]
		
	ops.beamIntegration('Lobatto', 16, secTag, 2)
	Nodes=[5000,20000+(numX-2)/2]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)


	CT=101010102
	eleNo=CT
	beamNormal=[1,0,0]
	Nodes=[20000+(numX-2)/2,6000]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+1+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+1+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+1+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+1+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)




	CT=101010105
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[6000,21000+(numX-2)/2]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+3+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+3+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+3+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+3+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)



	CT=101010106
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[21000+(numX-2)/2,7000]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+2+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+2+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+2+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+2+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)

	



	CT=101010103
	eleNo=CT
	beamNormal=[1,0,0]
	Nodes=[7000,22000+(numX-2)/2]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+4+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+4+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+4+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+4+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)




	CT=101010104
	eleNo=CT
	beamNormal=[1,0,0]
	Nodes=[22000+(numX-2)/2,8000]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+5+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+5+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+5+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+5+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)



	CT=101010107
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[8000,23000+(numX-2)/2]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+6+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+6+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+6+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+6+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)




	CT=101010108
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[23000+(numX-2)/2,5000]
	ops.geomTransf(coordTransf, CT, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', pc*1000+7+nodeNum, Nodes[0], Nodes[1], A, E, G, J, Iy, Iz, nodeNum+100000)
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', pc*1000+7+nodeNum, *[Nodes[0], Nodes[1]],CT, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+7+nodeNum, *[Nodes[0], Nodes[1]], CT, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', pc*1000+7+nodeNum, *[Nodes[0], Nodes[1]],3,secTag,CT,'-iter', 1)





if useSlabs==1:
	#  SLAB CORNER CONNECTIONS

	SN=np.array(np.linspace(100000,100000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))



	SN=np.reshape(SN,((numX+1),(numY+1)))


	LS1=SN[0][0]
	print(LS1)
	LS2=SN[0][numX]
	print(LS2)
	LS3=SN[numX][numX]
	print(LS3)
	LS4=SN[numX][0]
	print(LS4)

	FOAMySeesInstance.coupledNodes.remove(LS1)
	FOAMySeesInstance.coupledNodes.remove(LS2)
	FOAMySeesInstance.coupledNodes.remove(LS3)
	FOAMySeesInstance.coupledNodes.remove(LS4)
	eqDOFS=[[nElem,LS1],[1000+nElem,LS2],[2000+nElem,LS3],[3000+nElem,LS4]]
	for eqDOFpair in eqDOFS:
		rNodeTag=eqDOFpair[0]
		cNodeTag=eqDOFpair[1]
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		
		#ops.rigidLink('bar', rNodeTag, cNodeTag)
		#ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)


	SN=np.array(np.linspace(200000,200000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))



	SN=np.reshape(SN,((numX+1),(numY+1)))


	LS1=SN[0][0]
	print(LS1)
	LS2=SN[0][numX]
	print(LS2)
	LS3=SN[numX][numX]
	print(LS3)
	LS4=SN[numX][0]
	print(LS4)

	eqDOFS=[[5000+nElem2nd,LS1],[6000+nElem2nd,LS2],[7000+nElem2nd,LS3],[8000+nElem2nd,LS4]]
	for eqDOFpair in eqDOFS:
		rNodeTag=eqDOFpair[0]
		cNodeTag=eqDOFpair[1]
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
		
		#ops.rigidLink('bar', rNodeTag, cNodeTag)
		#ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)


	SN=np.array(np.linspace(300000,300000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))


	SN=np.reshape(SN,((numX+1),(numY+1)))


	LS1=SN[0][0]
	print(LS1)
	LS2=SN[0][numX]
	print(LS2)
	LS3=SN[numX][numX]
	print(LS3)
	LS4=SN[numX][0]
	print(LS4)

	eqDOFS=[[15000+nElemTOP,LS1],[16000+nElemTOP,LS2],[17000+nElemTOP,LS3],[18000+nElemTOP,LS4]]

	for eqDOFpair in eqDOFS:
		rNodeTag=eqDOFpair[0]
		cNodeTag=eqDOFpair[1]
		if constraintType=='EQDOF': 
			ops.equalDOF(rNodeTag, cNodeTag, *dofs)
		elif constraintType=='RIGIDLINK':	
			ops.rigidLink('bar', rNodeTag, cNodeTag)
		elif constraintType=='EQDOFMIXED':
			ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	
print(FOAMySeesInstance.coupledNodes)

if BaseIsolated==1:
	###############################################################################################################################################
	#   BASE ISOLATOR ELEMENTS #
	###############################################################################################################################################



	ops.node(9000000, *[Xmin,Ymin,0.0])
	ops.node(9001000, *[Xmin,Ymax,0.0])
	ops.node(9002000, *[Xmax,Ymax,0.0])
	ops.node(9003000, *[Xmax,Ymin,0.0])


	# eleTag=			 #(int)	unique element object tag
	# eleNodes=		   #(list (int))	a list of two element nodes

	kInit=50000			  #(float)	initial elastic stiffness in local shear direction
	qd=500				 #(float)	characteristic strength
	alpha1=0.01			#(float)	post yield stiffness ratio of linear hardening component
	alpha2=0.1			 #(float)	post yield stiffness ratio of non-linear hardening component
	mu=2				#(float)	exponent of non-linear hardening component

	matTag=222
	E=100000
	ops.uniaxialMaterial('Elastic', matTag, E) #, eta=0.0, Eneg=E)

	PMatTag=222			#(int)	tag associated with previously-defined UniaxialMaterial in axial direction
	TMatTag=222			#(int)	tag associated with previously-defined UniaxialMaterial in torsional direction
	MyMatTag=222		   #(int)	tag associated with previously-defined UniaxialMaterial in moment direction around local y-axis
	MzMatTag=222		   #(int)	tag associated with previously-defined UniaxialMaterial in moment direction around local z-axis

	# x1 x2 x3 (float)	vector components in global coordinates defining local x-axis (optional)
	# y1 y2 y3 (float)	vector components in global coordinates defining local y-axis (optional)
	# sDratio (float)	shear distance from iNode as a fraction of the element length (optional, default = 0.5)
	# '-doRayleigh' (str)	to include Rayleigh damping from the bearing (optional, default = no Rayleigh damping contribution)
	# m (float)	element mass (optional, default = 0.0)



	ops.element('elastomericBearingPlasticity', 9000001, *[0, 9000000], kInit, qd, alpha1, alpha2, mu, '-P', PMatTag, '-T', TMatTag, '-My', MyMatTag, '-Mz', MzMatTag) #, <'-orient', <x1, x2, x3>, y1, y2, y3>, <'-shearDist', sDratio>, <'-doRayleigh'>, <'-mass', m>)
	ops.element('elastomericBearingPlasticity', 9000002, *[1000, 9001000], kInit, qd, alpha1, alpha2, mu, '-P', PMatTag, '-T', TMatTag, '-My', MyMatTag, '-Mz', MzMatTag) #, <'-orient', <x1, x2, x3>, y1, y2, y3>, <'-shearDist', sDratio>, <'-doRayleigh'>, <'-mass', m>)
	ops.element('elastomericBearingPlasticity', 9000003, *[2000, 9002000], kInit, qd, alpha1, alpha2, mu, '-P', PMatTag, '-T', TMatTag, '-My', MyMatTag, '-Mz', MzMatTag) #, <'-orient', <x1, x2, x3>, y1, y2, y3>, <'-shearDist', sDratio>, <'-doRayleigh'>, <'-mass', m>)
	ops.element('elastomericBearingPlasticity', 9000004, *[3000, 9003000], kInit, qd, alpha1, alpha2, mu, '-P', PMatTag, '-T', TMatTag, '-My', MyMatTag, '-Mz', MzMatTag) #, <'-orient', <x1, x2, x3>, y1, y2, y3>, <'-shearDist', sDratio>, <'-doRayleigh'>, <'-mass', m>)

	ops.fix(9000000,*[1,1,1,1,1,1])
	ops.fix(9001000,*[1,1,1,1,1,1])
	ops.fix(9002000,*[1,1,1,1,1,1])
	ops.fix(9003000,*[1,1,1,1,1,1])
else:
	pass

ops.fixZ(2.0,*[1,1,1,1,1,1])


if useSlabs==1:
	pass
	# SN=np.array(np.linspace(100000,100000+(numX+1)*(numY+1)-1,(numX+1)*(numY+1)))

	# rmList=[]

	# SN=np.reshape(SN,((numX+1),(numY+1)))
	# LS1=SN[0][1]

	# FOAMySeesInstance.coupledNodes.remove(LS1)
	# rmList.append(LS1)
	# LS1=SN[1][1]

	# FOAMySeesInstance.coupledNodes.remove(LS1)
	# rmList.append(LS1)
	# LS1=SN[1][0]

	# FOAMySeesInstance.coupledNodes.remove(LS1)
	# rmList.append(LS1)
	# LS2=SN[1][numX]

	# FOAMySeesInstance.coupledNodes.remove(LS2)
	# rmList.append(LS2)
	# LS2=SN[0][numX-1]
	
	# rmList.append(LS2)
	# FOAMySeesInstance.coupledNodes.remove(LS2)
	# LS2=SN[1][numX-1]
	# rmList.append(LS2)
	# FOAMySeesInstance.coupledNodes.remove(LS2)

	# LS3=SN[numX-1][numX-1]
	# rmList.append(LS3)
	# FOAMySeesInstance.coupledNodes.remove(LS3)

	# LS3=SN[numX][numX-1]
	# rmList.append(LS3)
	# FOAMySeesInstance.coupledNodes.remove(LS3)

	# LS3=SN[numX-1][numX]
	# rmList.append(LS3)
	# FOAMySeesInstance.coupledNodes.remove(LS3)


	# LS4=SN[numX][1]
	# rmList.append(LS4)
	# FOAMySeesInstance.coupledNodes.remove(LS4)

	# LS4=SN[numX-1][1]
	# rmList.append(LS4)
	# FOAMySeesInstance.coupledNodes.remove(LS4)

	# LS4=SN[numX-1][0]
	# rmList.append(LS4)
	# FOAMySeesInstance.coupledNodes.remove(LS4)


ops.recorder('Node', '-file', 'reactionNode0.out','-time', '-node', 0, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'reactionNode1000.out','-time', '-node', 1000, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'reactionNode2000.out','-time', '-node', 2000, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'reactionNode3000.out','-time', '-node', 3000, '-dof', 1,2,3,4,5,6, 'reaction')


numEigenvalues=1
print('Eigen ',ops.eigen('-genBandArpack', numEigenvalues))

# opsplt.createODB("3DFrame", "Gravity", Nmodes=3)

############################
# Gravity Analysis
############################

# print("running gravity")

# timeSeries("Linear", 1)					# create TimeSeries for gravity analysis
# pattern('Plain',1,1)
# load(10000,0,0,-200.0,0,0,0)	# apply vertical load

# # recorder('Node','-file','Disp.txt','-time','-node',5000,'-dof',1,'disp')

# constraints('Transformation')
# numberer('ParallelRCM')
# #numberer('Plain')
# #system('BandGen')
# system('Mumps')
# #test('NormDispIncr',1.0e-11,200)
# algorithm('Linear') 
# integrator('CentralDifference')
# analysis('Transient')
# #domainChange()
# print("gravity analysis started...")
# partition()   
# domainChange()
# for j in range(10):
	
#	 analyze(1,1e-10)
#	 domainChange()
# print("gravity analysis complete...")
# createODB("3DFrame", "Gravity")

# # plot_deformedshape(Model="3DFrame", LoadCase="Gravity")
# loadConst('-time',0.0)					# Keep the gravity loads for further analysis

# wipeAnalysis()
# AnalysisType="Pushover"
# ###############################
# ### Cyclic ANALYSIS
# ###############################

# if(AnalysisType=="Cyclic"):
	
#	 # This is a load controlled analysis. The input load file "RCshearwall_Load_input.txt" should be in the 
#	 # .. same folder as the model file.
	
#	 print("<<<< Running Cyclic Analysis >>>>")
	
#	 timeSeries('Path',2,'-dt',0.1,'-filePath','RCshearwall_Load_input.txt')
#	 pattern('Plain',2,2)
#	 sp(10000,1,1)								# construct a single-point constraint object added to the LoadPattern.

#	 constraints('Penalty',1e20,1e20)
#	 numberer('RCM')
#	 system('BandGeneral')
#	 test('NormDispIncr',1e-05, 100, 1)
#	 algorithm('KrylovNewton')
#	 integrator('LoadControl',0.1)
#	 analysis('Static')
#	 analyze(700)


# #######################
# # PUSHOVER ANALYSIS
# #######################

# if(AnalysisType=="Pushover"):
	
#	 print("<<<< Running Pushover Analysis >>>>")

#	 # create a plain load pattern for pushover analysis
#	 pattern("Plain", 2, 1)
	
#	 ControlNode=10000
#	 ControlDOF=1
#	 MaxDisp= 0.020
#	 DispIncr=0.00001
#	 NstepsPush=int(MaxDisp/DispIncr)
	
#	 load(ControlNode, 1.00, 0.0, 0.0,0,0,0)	# Apply a unit reference load in DOF=1
	
#	 system("BandGeneral")
#	 numberer("RCM")
#	 constraints('Penalty',1e20,1e20)
#	 integrator("DisplacementControl", ControlNode, ControlDOF, DispIncr)
#	 algorithm('KrylovNewton')
#	 test('NormDispIncr',1e-05, 1000, 2)
#	 analysis("Static")
	
#	 numEigenvalues=3
	

#	 # Create a folder to put the output
#	 PushDataDir = r'PushoverOut'
#	 if not os.path.exists(PushDataDir):
#		 os.makedirs(PushDataDir)
#	 recorder('Node', '-file', "PushoverOut/React.out", '-closeOnWrite', '-node', 1, 2, 3, 4, 5, '-dof',1, 'reaction')
#	 recorder('Node', '-file', "PushoverOut/Disp.out", '-closeOnWrite', '-node', ControlNode, '-dof',1, 'disp')

#	 # Perform pushover analysis
#	 dataPush = np.zeros((NstepsPush+1,5))
#	 for j in range(NstepsPush):
#		 analyze(1)
#		 dataPush[j+1,0] = nodeDisp(ControlNode,1)*1000		# Convert to mm
#		 dataPush[j+1,1] = -getLoadFactor(2)*0.001			# Convert to kN
		
#	 # Read test output data to plot
#	 Test = np.loadtxt("RCshearwall_TestOutput.txt", delimiter="\t", unpack="False")
	
#	 ## Set parameters for the plot
#	 plt.rcParams.update({'font.size': 7})
#	 plt.figure(figsize=(4,3), dpi=100)
#	 plt.rc('font', family='serif')
#	 plt.plot(Test[0,:], Test[1,:], color="black", linewidth=0.8, linestyle="--", label='Test')
#	 plt.plot(dataPush[:,0], -dataPush[:,1], color="red", linewidth=1.2, linestyle="-", label='Pushover')
#	 plt.axhline(0, color='black', linewidth=0.4)
#	 plt.axvline(0, color='black', linewidth=0.4)
#	 plt.xlim(-25, 25)
#	 plt.xticks(np.linspace(-20,20,11,endpoint=True)) 
#	 plt.grid(linestyle='dotted') 
#	 plt.xlabel('Displacement (mm)')
#	 plt.ylabel('Base Shear (kN)')
#	 plt.legend()
#	 plt.savefig("PushoverOut/RCshearwall_PushoverCurve.png",dpi=1200)
#	 plt.show()
	
	
	