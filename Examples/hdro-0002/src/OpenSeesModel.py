
ops.wipe()

##########################################


"""
Created on Tue Sep 26 16:15:31 2023

@author: nsaoirse
"""



##########################################################
#                                                         #
# Procedure to compute ultimate lateral resistance, p_u,  #
#  and displacement at 50% of lateral capacity, y50, for  #
#  p-y springs representing cohesionless soil.            #
#   Converted to openseespy by: Pavan Chigullapally       #
#                               University of Auckland    # 
#                                                         #
#   Created by:   Hyung-suk Shin                          #
#                 University of Washington                #
#   Modified by:  Chris McGann                            #
#                 Pedro Arduino                           #
#                 Peter Mackenzie-Helnwein                #
#                 University of Washington                #
#                                                         #
###########################################################

# references
#  American Petroleum Institute (API) (1987). Recommended Practice for Planning, Designing and
#   Constructing Fixed Offshore Platforms. API Recommended Practice 2A(RP-2A), Washington D.C,
#   17th edition.
#
# Brinch Hansen, J. (1961). "The ultimate resistance of rigid piles against transversal forces."
#  Bulletin No. 12, Geoteknisk Institute, Copenhagen, 59.
#
#  Boulanger, R. W., Kutter, B. L., Brandenberg, S. J., Singh, P., and Chang, D. (2003). Pile 
#   Foundations in liquefied and laterally spreading ground during earthquakes: Centrifuge experiments
#   and analyses. Center for Geotechnical Modeling, University of California at Davis, Davis, CA.
#   Rep. UCD/CGM-03/01.
#
#  Reese, L.C. and Van Impe, W.F. (2001), Single Piles and Pile Groups Under Lateral Loading.
#    A.A. Balkema, Rotterdam, Netherlands.

import math

def get_pyParam( pyDepth, gamma, phiDegree, b, pEleLength, puSwitch, kSwitch, gwtSwitch):
    
    #----------------------------------------------------------
    #  define ultimate lateral resistance, pult 
    #----------------------------------------------------------
    
    # pult is defined per API recommendations (Reese and Van Impe, 2001 or API, 1987) for puSwitch = 1
    #  OR per the method of Brinch Hansen (1961) for puSwitch = 2
    
    pi = 3.14159265358979
    phi = phiDegree * (pi/180)
    zbRatio = pyDepth / b
    
    #-------API recommended method-------
    
    if puSwitch == 1:
    
      # obtain loading-type coefficient A for given depth-to-diameter ratio zb
      #  ---> values are obtained from a figure and are therefore approximate
        zb = []
        dataNum = 41
        for i in range(dataNum):
            b1 = i * 0.125
            zb.append(b1)
        As = [2.8460, 2.7105, 2.6242, 2.5257, 2.4271, 2.3409, 2.2546, 2.1437, 2.0575, 1.9589, 1.8973, 1.8111, 1.7372, 1.6632, 1.5893, 1.5277, 1.4415, 1.3799, 1.3368, 1.2690, 1.2074, 1.1581, 
            1.1211, 1.0780, 1.0349, 1.0164, 0.9979, 0.9733, 0.9610, 0.9487, 0.9363, 0.9117, 0.8994, 0.8994, 0.8871, 0.8871, 0.8809, 0.8809, 0.8809, 0.8809, 0.8809] 
      
      # linear interpolation to define A for intermediate values of depth:diameter ratio
        for i in range(dataNum):
            if zbRatio >= 5.0:
                A = 0.88
            elif zb[i] <= zbRatio and zbRatio <= zb[i+1]:
                A = (As[i+1] - As[i])/(zb[i+1] - zb[i]) * (zbRatio-zb[i]) + As[i]
                
      # define common terms
        alpha = phi / 2
        beta = pi / 4 + phi / 2
        K0 = 0.4
        
        tan_1 = math.tan(pi / 4 - phi / 2)        
        Ka = math.pow(tan_1 , 2) 
    
      # terms for Equation (3.44), Reese and Van Impe (2001)
        tan_2 = math.tan(phi)
        tan_3 = math.tan(beta - phi)
        sin_1 = math.sin(beta)
        cos_1 = math.cos(alpha)
        c1 = K0 * tan_2 * sin_1 / (tan_3*cos_1)
        
        tan_4 = math.tan(beta)
        tan_5 = math.tan(alpha)
        c2 = (tan_4/tan_3)*tan_4 * tan_5
        
        c3 = K0 * tan_4 * (tan_2 * sin_1 - tan_5)
        
        c4 = tan_4 / tan_3 - Ka
    
        # terms for Equation (3.45), Reese and Van Impe (2001)
        pow_1 = math.pow(tan_4,8)
        pow_2 = math.pow(tan_4,4)
        c5 = Ka * (pow_1-1)
        c6 = K0 * tan_2 * pow_2
    
      # Equation (3.44), Reese and Van Impe (2001)
        pst = gamma * pyDepth * (pyDepth * (c1 + c2 + c3) + b * c4)
    
      # Equation (3.45), Reese and Van Impe (2001)
        psd = b * gamma * pyDepth * (c5 + c6)
    
      # pult is the lesser of pst and psd. At surface, an arbitrary value is defined
        if pst <=psd:
            if pyDepth == 0:
                pu = 0.01
              
            else:
                pu = A * pst
              
        else:
            pu = A * psd
          
      # PySimple1 material formulated with pult as a force, not force/length, multiply by trib. length
        pult = pu * pEleLength
    
    #-------Brinch Hansen method-------
    elif puSwitch == 2:
      # pressure at ground surface
        cos_2 = math.cos(phi)
        
        tan_6 = math.tan(pi/4+phi/2) 
        
        sin_2 = math.sin(phi)
        sin_3 = math.sin(pi/4 + phi/2)
        
        exp_1 = math.exp((pi/2+phi)*tan_2)
        exp_2 = math.exp(-(pi/2-phi) * tan_2)
        
        Kqo = exp_1 * cos_2 * tan_6 - exp_2 * cos_2 * tan_1
        Kco = (1/tan_2) * (exp_1 * cos_2 * tan_6 - 1)
    
      # pressure at great depth
        exp_3 = math.exp(pi * tan_2)
        pow_3 = math.pow(tan_2,4)
        pow_4 = math.pow(tan_6,2)
        dcinf = 1.58 + 4.09 * (pow_3)
        Nc = (1/tan_2)*(exp_3)*(pow_4 - 1)
        Ko = 1 - sin_2
        Kcinf = Nc * dcinf
        Kqinf = Kcinf * Ko * tan_2
    
      # pressure at an arbitrary depth
        aq = (Kqo/(Kqinf - Kqo))*(Ko*sin_2/sin_3)
        KqD = (Kqo + Kqinf * aq * zbRatio)/(1 + aq * zbRatio)
    
      # ultimate lateral resistance
        if pyDepth == 0:
            pu = 0.01
        else:
            pu = gamma * pyDepth * KqD * b
               
      # PySimple1 material formulated with pult as a force, not force/length, multiply by trib. length
        pult  = pu * pEleLength
        
    #----------------------------------------------------------
    #  define displacement at 50% lateral capacity, y50
    #----------------------------------------------------------
    
    # values of y50 depend of the coefficent of subgrade reaction, k, which can be defined in several ways.
    #  for gwtSwitch = 1, k reflects soil above the groundwater table
    #  for gwtSwitch = 2, k reflects soil below the groundwater table
    #  a linear variation of k with depth is defined for kSwitch = 1 after API (1987)
    #  a parabolic variation of k with depth is defined for kSwitch = 2 after Boulanger et al. (2003)
    
    # API (1987) recommended subgrade modulus for given friction angle, values obtained from figure (approximate)
    
    ph = [28.8, 29.5, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0]    
   
    # subgrade modulus above the water table
    if gwtSwitch == 1:
        k = [10, 23, 45, 61, 80, 100, 120, 140, 160, 182, 215, 250, 275]
        
    else:
        k = [10, 20, 33, 42, 50, 60, 70, 85, 95, 107, 122, 141, 155]
    
    dataNum = 13  
    for i in range(dataNum):
        if ph[i] <= phiDegree and phiDegree <= ph[i+1]:
            khat = (k[i+1]-k[i])/(ph[i+1]-ph[i])*(phiDegree - ph[i]) + k[i]            
            
    # change units from (lb/in^3) to (kN/m^3)
    k_SIunits = khat * 271.45
    
    # define parabolic distribution of k with depth if desired (i.e. lin_par switch == 2)
    sigV = pyDepth * gamma
    
    if sigV == 0:
         sigV = 0.01
         
    if kSwitch == 2:
       # Equation (5-16), Boulanger et al. (2003)
        cSigma = math.pow(50 / sigV , 0.5)
       # Equation (5-15), Boulanger et al. (2003)
        k_SIunits = cSigma * k_SIunits
    
    # define y50 based on pult and subgrade modulus k
    
    # based on API (1987) recommendations, p-y curves are described using tanh functions.
    #  tcl does not have the atanh function, so must define this specifically
    
    #  i.e.  atanh(x) = 1/2*ln((1+x)/(1-x)), |x| < 1
    
    # when half of full resistance has been mobilized, p(y50)/pult = 0.5
    x = 0.5
    log_1 = math.log((1+x)/(1-x))
    atanh_value = 0.5 * log_1
    
    # need to be careful at ground surface (don't want to divide by zero)
    if pyDepth == 0.0:
        pyDepth = 0.01

    y50 = 0.5 * (pu/ A)/(k_SIunits * pyDepth) * atanh_value
    # return pult and y50 parameters
    outResult = []
    outResult.append(pult)
    outResult.append(y50)
    
    return outResult

#########################################################################################################################################################################

#########################################################################################################################################################################

###########################################################
#                                                         #
# Procedure to compute ultimate tip resistance, qult, and #
#  displacement at 50% mobilization of qult, z50, for     #
#  use in q-z curves for cohesionless soil.               #
#   Converted to openseespy by: Pavan Chigullapally       #  
#                               University of Auckland    #
#   Created by:  Chris McGann                             #
#                Pedro Arduino                            #
#                University of Washington                 #
#                                                         #
###########################################################

# references
#  Meyerhof G.G. (1976). "Bearing capacity and settlement of pile foundations." 
#   J. Geotech. Eng. Div., ASCE, 102(3), 195-228.
#
#  Vijayvergiya, V.N. (1977). "Load-movement characteristics of piles."
#   Proc., Ports 77 Conf., ASCE, New York.
#
#  Kulhawy, F.H. ad Mayne, P.W. (1990). Manual on Estimating Soil Properties for 
#   Foundation Design. Electrical Power Research Institute. EPRI EL-6800, 
#   Project 1493-6 Final Report.

def get_qzParam(phiDegree, b, sigV, G):
    
    # define required constants; pi, atmospheric pressure (kPa), pa, and coeff. of lat earth pressure, Ko
    pi = 3.14159265358979
    pa = 101
    sin_4 = math.sin(phiDegree * (pi/180))
    Ko = 1 - sin_4

  # ultimate tip pressure can be computed by qult = Nq*sigV after Meyerhof (1976)
  #  where Nq is a bearing capacity factor, phi is friction angle, and sigV is eff. overburden
  #  stress at the pile tip.
    phi = phiDegree * (pi/180)

  # rigidity index
    tan_7 = math.tan(phi)
    Ir = G/(sigV * tan_7)
  # bearing capacity factor
    tan_8 = math.tan(pi/4+phi/2)
    sin_5 = math.sin(phi)
    pow_4 = math.pow(tan_8,2)
    pow_5 = math.pow(Ir,(4*sin_5)/(3*(1+sin_5)))
    exp_4 = math.exp(pi/2-phi)
    
    Nq = (1+2*Ko)*(1/(3-sin_5))*exp_4*(pow_4)*(pow_5)  
  # tip resistance
    qu = Nq * sigV
  # QzSimple1 material formulated with qult as force, not stress, multiply by area of pile tip
    pow_6 = math.pow(b, 2)  
    qult = qu * pi*pow_6/4

  # the q-z curve of Vijayvergiya (1977) has the form, q(z) = qult*(z/zc)^(1/3)
  #  where zc is critical tip deflection given as ranging from 3-9% of the
  #  pile diameter at the tip.  

  # assume zc is 5% of pile diameter
    zc = 0.05 * b

  # based on Vijayvergiya (1977) curve, z50 = 0.125*zc
    z50 = 0.125 * zc

  # return values of qult and z50 for use in q-z material
    outResult = []
    outResult.append(qult)
    outResult.append(z50)
    
    return outResult

#########################################################################################################################################################################

#########################################################################################################################################################################
##########################################################
#                                                         #
# Procedure to compute ultimate resistance, tult, and     #
#  displacement at 50% mobilization of tult, z50, for     #
#  use in t-z curves for cohesionless soil.               #
#   Converted to openseespy by: Pavan Chigullapally       #
#                               University of Auckland    #
#   Created by:  Chris McGann                             #
#                University of Washington                 #
#                                                         #
###########################################################

def get_tzParam( phi, b, sigV, pEleLength):

# references
#  Mosher, R.L. (1984). "Load transfer criteria for numerical analysis of
#   axial loaded piles in sand." U.S. Army Engineering and Waterways
#   Experimental Station, Automatic Data Processing Center, Vicksburg, Miss.
#
#  Kulhawy, F.H. (1991). "Drilled shaft foundations." Foundation engineering
#   handbook, 2nd Ed., Chap 14, H.-Y. Fang ed., Van Nostrand Reinhold, New York

    pi = 3.14159265358979
    
  # Compute tult based on tult = Ko*sigV*pi*dia*tan(delta), where
  #   Ko    is coeff. of lateral earth pressure at rest, 
  #         taken as Ko = 0.4
  #   delta is interface friction between soil and pile,
  #         taken as delta = 0.8*phi to be representative of a 
  #         smooth precast concrete pile after Kulhawy (1991)
  
    delta = 0.8 * phi * pi/180

  # if z = 0 (ground surface) need to specify a small non-zero value of sigV
  
    if sigV == 0.0:
        sigV = 0.01
    
    tan_9 = math.tan(delta)
    tu = 0.4 * sigV * pi * b * tan_9
    
  # TzSimple1 material formulated with tult as force, not stress, multiply by tributary length of pile
    tult = tu * pEleLength

  # Mosher (1984) provides recommended initial tangents based on friction angle
	# values are in units of psf/in
    kf = [6000, 10000, 10000, 14000, 14000, 18000]
    fric = [28, 31, 32, 34, 35, 38]

    dataNum = len(fric)
    
    
	# determine kf for input value of phi, linear interpolation for intermediate values
    if phi < fric[0]:
        k = kf[0]
    elif phi > fric[5]:
        k = kf[5]
    else:
        for i in range(dataNum):
            if fric[i] <= phi and phi <= fric[i+1]:
                k = ((kf[i+1] - kf[i])/(fric[i+1] - fric[i])) * (phi - fric[i]) + kf[i]
        

  # need to convert kf to units of kN/m^3
    kSIunits =  k * 1.885

  # based on a t-z curve of the shape recommended by Mosher (1984), z50 = tult/kf
    z50 = tult / kSIunits

  # return values of tult and z50 for use in t-z material
    outResult = []
    outResult.append(tult)
    outResult.append(z50)

    return outResult


#########################################################################################################################################################################

#########################################################################################################################################################################

###########################################################
#                                                         #
# Static pushover of a single pile, modeled as a beam on  #
#  a nonlinear Winkler foundation.  Lateral soil response #
#  is described by p-y springs.  Vertical soil response   #
#  described by t-z and q-z springs.                      #
#   Converted to openseespy by: Pavan Chigullapally       #
#                               University of Auckland    #
#   Created by:  Chris McGann                             #
#                HyungSuk Shin                            #
#                Pedro Arduino                            #
#                Peter Mackenzie-Helnwein                 #
#              --University of Washington--               #
#                                                         #
# ---> Basic units are kN and meters                      #
#                                                         #
###########################################################

ops.wipe()

#########################################################################################################################################################################

#########################################################################################################################################################################


#########################################


FOAMySeesInstance.coupledNodes=[]

nElemBOTTOM=16
nElem2nd=4
nElemTOP=4

numX=2
numY=2

ColumnConnectionLength=0.05
BeamConnLength=0.0254*4

BuildingLocation=[(40.88230+41.89830)/2, 0.0]
BuildingDepth=41.89830-40.88230
BuildingWidth=0.508*2
Xmin=BuildingLocation[0]-BuildingDepth/2
Xmax=BuildingLocation[0]+BuildingDepth/2
Ymin=BuildingLocation[1]-BuildingWidth/2
Ymax=BuildingLocation[1]+BuildingWidth/2

zPileBase=1.0
zBase=2.0
z1stFloor=2.6180000
z2ndFloor=2.6180000+(27.1*0.0254)
z3rdFloor=2.6180000+(27.1*0.0254)+(18*0.0254)
z4thFloor=2.6180000+(27.1*0.0254)+(18*0.0254)*2

beamType='NLBeamCol'
beamType='forceBeamCol'


beamType='elastic'
beamType='dispBeamCol'






noBraces=0
useBeams=1

constraintType='EQDOFMIXED'
constraintType='RIGIDLINK'	
constraintType='EQDOF'

useSlabs=1

BaseIsolated=0

zBumpDebug=0.0

numDOF=6
rcdofs=[1,1,2,2,3,3,4,4,5,5,6,6]
dofs=[1,2,3,4,5,6]
soildofs=[1,3]
FOAMySeesInstance.osi=ops.model('basic','-ndm',3,'-ndf',numDOF)
# all the units are in SI units N and mm

#----------------------------------------------------------
#  soil properties
#----------------------------------------------------------

# soil unit weight (kN/m^3)
gamma = 17.0
# soil internal friction angle (degrees)
phi = 36.0
# soil shear modulus at pile tip (kPa)
Gsoil = 150000.0

# select pult definition method for p-y curves
# API (default) --> 1
# Brinch Hansen --> 2
puSwitch = 1 

# variation in coefficent of subgrade reaction with depth for p-y curves
# API linear variation (default)   --> 1
# modified API parabolic variation --> 2
kSwitch = 1

# effect of ground water on subgrade reaction modulus for p-y curves
# above gwt --> 1
# below gwt --> 2
gwtSwitch = 1
#----------------------------------------------------------
#  pile geometry and mesh
#----------------------------------------------------------

# length of pile head (above ground surface) (m)
L1 = 0.0
# length of embedded pile (below ground surface) (m)
L2 = 1.0
# pile diameter
diameter = 4*0.0254

# number of pile elements
nElePile = nElemBOTTOM
# pile element length 
eleSize = (L1+L2)/nElePile

# number of total pile nodes
nNodePile =  1 + nElePile

  # vertical effective stress at pile tip, no water table (depth is embedded pile length)
sigVq = gamma * L2
  # procedure to define qult and z50
qzParam = get_qzParam(phi, diameter, sigVq, Gsoil)
qult = qzParam [0]
z50q = qzParam [1]

#ops.uniaxialMaterial('QzSimple1', 101, 2, qult, z50q) #, 0.0, 0.0
ops.uniaxialMaterial('TzSimple1', 101, 2, qult, z50q, 0.0)

# ------------------------------
# Start of analysis generation
# ------------------------------


A=0.49*0.00064516
Iz=0.118*0.00064516*0.00064516
Iy=0.118*0.00064516*0.00064516
Jxx=0.237*0.00064516*0.00064516
E=29000*6895000 #ksi*conversion
G=E/(2*(1.3))
secTag=101
ops.section('Elastic', secTag, E, A, Iz, Iy, G, Jxx)

ops.beamIntegration('Lobatto', 3,secTag, 2)
# section('Elastic', secTag, E_mod, A, Iz, G_mod=None)
A=0.84*0.00064516
Iz=0.486*0.00064516*0.00064516
Iy=0.486*0.00064516*0.00064516
Jxx=0.796*0.00064516*0.00064516
E=29000*6895000 #ksi*conversion
G=E/(2*(1.3))
secTag=15
ops.section('Elastic', secTag, E, A, Iz, Iy, G, Jxx)


ops.beamIntegration('Lobatto', 2, secTag, 2)

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
OuterRad=2
OuterDiam=OuterRad*2
CFTWallT=0.25

## Sections]
secTag=501
ops.section('Fiber', secTag, '-GJ', 10000000000.00000)

matTag=650
numSubdivCirc=8
numSubdivRad=8
center=[0.00000, 0.00000]
rad=[0.00000, (OuterRad-CFTWallT)*0.0254] #INCH*conversionToMeters
ang=[0.00000, 360.00000]
ops.patch('circ', matTag, numSubdivCirc, numSubdivRad, *center, *rad, *ang)

matTag=600
numSubdivCirc=8
numSubdivRad=4
center=[0.00000, 0.00000]
rad=[(OuterRad-CFTWallT)*0.0254, OuterRad*0.0254] #INCH*conversionToMeters
ang=[0.00000, 360.00000]
ops.patch('circ', matTag, numSubdivCirc, numSubdivRad, *center, *rad, *ang)
ops.beamIntegration('Lobatto', 1, secTag, 2)




# below Grade Beams

# SUB STORY COLUMNS
##################################################################################################################################################################


slabMass=1000

nElem=nElemBOTTOM
Zbase=zPileBase+ColumnConnectionLength
Ztop=zBase-ColumnConnectionLength
sectionLength=Ztop-Zbase
pileMass=sectionLength*7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*3.141519

subStoryColumnNodes=[]

subStoryColumnBaseNodes=[]
pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]

ops.node(123456780,Xmin,Ymin,1.0) 
ops.node(123456781,Xmin,Ymax,1.0)
ops.node(123456782,Xmax,Ymax,1.0)
ops.node(123456783,Xmax,Ymin,1.0)
bottombottomnodes=[123456780,123456781,123456782,123456783]
pc=100
startElem=90000000
startNode=90000000
soilnodes=[]
for location in pileLocs:

	node1=[location[0], location[1],  Zbase]
	node2=[location[0], location[1], Ztop]

	beamNormal=[-1,0,0]
	startElem+=pc
	startNode+=pc
	count=0
	xNodeList=np.linspace(node1[0],node2[0],nElemBOTTOM+1)
	yNodeList=np.linspace(node1[1],node2[1],nElemBOTTOM+1)
	zNodeList=np.linspace(node1[2],node2[2],nElemBOTTOM+1)
	nodalMass=(pileMass)/len(xNodeList)

	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.node(nodeNum, xNodeList[nodeNum-pc*1000],yNodeList[nodeNum-pc*1000],zNodeList[nodeNum-pc*1000])
		
        
		# ELENo=startElem
		# NN1=startNode+count
		
		# NN2=startNode+count+1
		# count+=2
        # #----------------------------------------------------------
        # #  create spring material objects
        # #----------------------------------------------------------
        
           
        # # t-z spring material    
		# pyDepth=zBase-zNodeList[nodeNum-pc*1000]
            # # vertical effective stress at current depth    
		# sigV = gamma * pyDepth
            # # procedure to define tult and z50
		# tzParam = get_tzParam(phi, diameter, sigV, eleSize)
		# tult = tzParam [0]
		# z50 = tzParam [1]
		# ops.uniaxialMaterial('TzSimple1', NN2, 2, tult, z50, 0.0)
        # # p-y spring material
        # # procedure to define pult and y50
		# pyParam = get_pyParam(pyDepth, gamma, phi, diameter, eleSize, puSwitch, kSwitch, gwtSwitch)
		# pult = pyParam [0]
		# y50 = pyParam [1]    
		# ops.uniaxialMaterial('PySimple1', NN1, 2, pult, y50, 0.0)       
        
		# ops.node(NN1,location[0], location[1], pyDepth)
		# ops.node(NN2,location[0], location[1], pyDepth)
		# ops.fix(NN1, 1, 1, 1, 1, 1, 1)
		# ops.fix(NN2, 0, 0, 1, 1, 1, 1)
		# ops.element('zeroLength', NN1, NN1, NN2, '-mat', NN1, NN2, '-dir', 1, 3)
		# ops.element('zeroLength', NN2, NN1, NN2, '-mat', NN1, NN2, '-dir', 2, 3)
		# ops.equalDOF(nodeNum, NN2, *dofs)
		# soilnodes.append(NN1)
		# soilnodes.append(NN2)
        

        
        
        #FOAMySeesInstance.coupledNodes.append(nodeNum)
	coordTransf = 'Corotational'

	coordTransf='Linear'
	coordTransf='PDelta'
	#############################

	## Model
	subStoryColumnBaseNodes.append(pc*1000)

	subStoryColumnNodes.append(pc*1000+len(xNodeList)-1)

	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, 501, pc*1000+nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1,'-mass', 3000)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1
ct=0




# FIRST STORY COLUMNS
##################################################################################################################################################################

slabMass=1000

nElem=nElemBOTTOM
Zbase=zBase
Ztop=z1stFloor
sectionLength=Ztop-Zbase
pileMass=sectionLength*7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*3.141519


firstStoryColumnNodes=[]
firstStoryColumnBaseNodes=[]

pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]




pc=1000
for location in pileLocs:

	node1=[location[0], location[1],  Zbase]
	node2=[location[0], location[1], Ztop]

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
	firstStoryColumnBaseNodes.append(pc*1000)

	firstStoryColumnNodes.append(pc*1000+len(xNodeList)-1)

	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, 501, pc*1000+nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1,'-mass', 3000)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1
    
ct=0

StrainGaugeElements=[]


bottombottomnodes=[123456780,123456781,123456782,123456783]

# beamEndNode=firstSlabCornerNodeNums[-1][0]
for corner in bottombottomnodes:
	columnTopNode=subStoryColumnBaseNodes[ct]
	print(columnTopNode,'columnTopNode')
	columnBottomNode=corner
	print(columnBottomNode,'columnBottomNode')
	ct+=1	
	beamNormal=[-1,0,0]
	ops.geomTransf(coordTransf, pc*1000+ct+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	nodeNum+=1

	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, columnBottomNode, columnTopNode, 501, pc*1000+ct+100000)  		
	elif beamType=='dispBeamCol':
  		ops.element('dispBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],pc*1000+ct+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
  		ops.element('forceBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode], pc*1000+ct+100000, 1)		
	elif beamType=='NLBeamCol':
  		ops.element('nonlinearBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],3,secTag,pc*1000+ct+100000,'-iter', 1)
	StrainGaugeElements.append(nodeNum)
ct=0
pc+=1
# beamEndNode=firstSlabCornerNodeNums[-1][0]
for corner in subStoryColumnNodes:
	columnTopNode=firstStoryColumnBaseNodes[ct]
	print(columnTopNode,'columnTopNode')
	columnBottomNode=corner
	print(columnBottomNode,'columnBottomNode')
	ct+=1	
	beamNormal=[-1,0,0]
	ops.geomTransf(coordTransf, pc*1000+ct+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	nodeNum+=1

	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, columnBottomNode, columnTopNode, 501, pc*1000+ct+100000)  		
	elif beamType=='dispBeamCol':
  		ops.element('dispBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],pc*1000+ct+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
  		ops.element('forceBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode], pc*1000+ct+100000, 1)		
	elif beamType=='NLBeamCol':
  		ops.element('nonlinearBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],3,secTag,pc*1000+ct+100000,'-iter', 1)
	StrainGaugeElements.append(nodeNum)

numEigenvalues=1
#print('Eigen ',ops.eigen('-genBandArpack', numEigenvalues))

# render the model after defining all the nodes and elements
#vfo.plot_model()




firstStoryBeamEndNodes=[]
	
	

print('I am here beams are made')

#vfo.plot_model()


# 2nd STORY COLUMNS
##################################################################################################################################################################



Zbase=z1stFloor

Ztop=z2ndFloor
sectionLength=Ztop-Zbase
pileMass=sectionLength*7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*3.141519


secondStoryColumnTopNodes=[]

secondStoryColumnBaseNodes=[]
pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]

pc=5
for location in pileLocs:

	node1=[location[0], location[1],  Zbase+ColumnConnectionLength]
	node2=[location[0], location[1], Ztop]

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
	secondStoryColumnBaseNodes.append(pc*1000)
	secondStoryColumnTopNodes.append(pc*1000+len(xNodeList)-1)



	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, 501, pc*1000+nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1

nodeNum+=1
ct=0
secondStoryColumnBaseNodes.append(secondStoryColumnBaseNodes[0])
secondStoryColumnTopNodes.append(secondStoryColumnTopNodes[0])

# beamEndNode=firstSlabCornerNodeNums[-1][0]
for corner in firstStoryColumnNodes:
	columnTopNode=secondStoryColumnBaseNodes[ct]
	print(columnTopNode,'columnTopNode')
	columnBottomNode=corner
	print(columnBottomNode,'columnBottomNode')
	ct+=1	
	beamNormal=[-1,0,0]
	ops.geomTransf(coordTransf, pc*1000+ct+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	nodeNum+=1

	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, columnBottomNode, columnTopNode, 501, pc*1000+ct+100000)  		
	elif beamType=='dispBeamCol':
  		ops.element('dispBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],pc*1000+ct+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
  		ops.element('forceBeamColumn', pc*1000+ct, *[ columnBottomNode, columnTopNode], pc*1000+ct+100000, 1)		
	elif beamType=='NLBeamCol':
  		ops.element('nonlinearBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],3,secTag,pc*1000+ct+100000,'-iter', 1)





ct+=1


numEigenvalues=1
# print('Eigen ',ops.eigen('-genBandArpack', numEigenvalues))

# render the model after defining all the nodes and elements
#vfo.plot_model()






integTag=15
	


thirdStoryColumnTopNodes=[]

thirdStoryColumnBaseNodes=[]
	
ops.beamIntegration('Lobatto', integTag, secTag, 2)
slabMass=1000

# 3rd STORY COLUMNS
##################################################################################################################################################################


Zbase=z2ndFloor+ColumnConnectionLength
#Ztop=1.8288
Ztop=z3rdFloor
sectionLength=Ztop-Zbase
pileMass=sectionLength*7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*3.141519


thirdStoryColumnTopNodes=[]

thirdStoryColumnBaseNodes=[]
pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]

pc=25
for location in pileLocs:

	node1=[location[0], location[1],  Zbase]
	node2=[location[0], location[1], Ztop]

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
	thirdStoryColumnBaseNodes.append(pc*1000)
	thirdStoryColumnTopNodes.append(pc*1000+len(xNodeList)-1)



	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, 501, pc*1000+nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1

nodeNum+=1
ct=0
thirdStoryColumnBaseNodes.append(thirdStoryColumnBaseNodes[0])
thirdStoryColumnTopNodes.append(thirdStoryColumnTopNodes[0])

# beamEndNode=firstSlabCornerNodeNums[-1][0]
for corner in secondStoryColumnTopNodes:
	columnTopNode=thirdStoryColumnBaseNodes[ct]
	print(columnTopNode,'columnTopNode')
	columnBottomNode=corner
	print(columnBottomNode,'columnBottomNode')
	ct+=1	
	beamNormal=[-1,0,0]
	ops.geomTransf(coordTransf, pc*1000+ct+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	nodeNum+=1

	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, columnBottomNode, columnTopNode, 501, pc*1000+ct+100000)  		
	elif beamType=='dispBeamCol':
  		ops.element('dispBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],pc*1000+ct+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
  		ops.element('forceBeamColumn', pc*1000+ct, *[ columnBottomNode, columnTopNode], pc*1000+ct+100000, 1)		
	elif beamType=='NLBeamCol':
  		ops.element('nonlinearBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],3,secTag,pc*1000+ct+100000,'-iter', 1)


#vfo.plot_model()



#vfo.plot_model()
# TOP STORY COLUMNS
##################################################################################################################################################################

Zbase=z3rdFloor+ColumnConnectionLength
#Ztop=1.8288
Ztop=z4thFloor
sectionLength=Ztop-Zbase
pileMass=sectionLength*7840*((0.5*OuterDiam*0.0254)**2 - (0.5*(OuterDiam-CFTWallT)*0.0254)**2)*(0.6180000)*3.141519 + 2400*((0.5*(OuterDiam-CFTWallT)*0.0254)**2)*3.141519



topStoryColumnBaseNodes=[]

topStoryColumnTopNodes=[]


pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]

pc=35
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
	topStoryColumnBaseNodes.append(pc*1000)
	topStoryColumnTopNodes.append(pc*1000+len(xNodeList)-1)



	nodRotMass=0.
	for nodeNum in range(pc*1000, pc*1000+len(xNodeList)):
		ops.mass(nodeNum,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])

	for nodeNum in range(1 + pc*1000, pc*1000+len(xNodeList)):
		ops.geomTransf(coordTransf, pc*1000+nodeNum+100000, beamNormal[0],beamNormal[1],beamNormal[2])
		
		if beamType=='elastic':
			ops.element('elasticBeamColumn', nodeNum, nodeNum-1, nodeNum, 501, pc*1000+nodeNum+100000)
		elif beamType=='dispBeamCol':
			ops.element('dispBeamColumn', nodeNum, *[nodeNum-1, nodeNum],pc*1000+nodeNum+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
		elif beamType=='forceBeamCol':
			ops.element('forceBeamColumn', pc*1000+nodeNum, *[nodeNum-1, nodeNum], pc*1000+nodeNum+100000, 1)		
		elif beamType=='NLBeamCol':
			ops.element('nonlinearBeamColumn', nodeNum, *[nodeNum-1, nodeNum],3,secTag,pc*1000+nodeNum+100000,'-iter', 1)
	pc+=1

ct=0
topStoryColumnBaseNodes.append(topStoryColumnBaseNodes[0])
topStoryColumnTopNodes.append(topStoryColumnTopNodes[0])

# beamEndNode=firstSlabCornerNodeNums[-1][0]
for corner in thirdStoryColumnTopNodes:
	columnTopNode=topStoryColumnBaseNodes[ct]
	print(columnTopNode,'columnTopNode')
	columnBottomNode=corner
	print(columnBottomNode,'columnBottomNode')
	ct+=1	
	beamNormal=[-1,0,0]
	ops.geomTransf(coordTransf, pc*1000+ct+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	nodeNum+=1

	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, columnBottomNode, columnTopNode, 501, pc*1000+ct+100000)  		
	elif beamType=='dispBeamCol':
  		ops.element('dispBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],pc*1000+ct+100000, 1)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
  		ops.element('forceBeamColumn', pc*1000+ct, *[ columnBottomNode, columnTopNode], pc*1000+ct+100000, 1)		
	elif beamType=='NLBeamCol':
  		ops.element('nonlinearBeamColumn', nodeNum, *[ columnBottomNode, columnTopNode],3,secTag,pc*1000+ct+100000,'-iter', 1)


BeamCenterNodesFirstFloor=[]
firstStoryColumnNodes.append(firstStoryColumnNodes[0])
ct=0
Gct=0
NnumMiddleStart=10000000
for nodeCurrent in firstStoryColumnNodes[0:4]:
	node2=firstStoryColumnNodes[ct+1]
	node1=nodeCurrent
	print(nodeCurrent)
	node1Loc=ops.nodeCoord(node1)
	node2Loc=ops.nodeCoord(node2)
    
    
    
	middleLoc=(np.array(node2Loc)+np.array(node1Loc))/2
    
	beamNormal=np.cross(np.array(node2Loc)-np.array(node1Loc),[0,0,1])
    

    
	ops.node(NnumMiddleStart, middleLoc[0],middleLoc[1],middleLoc[2])
	BeamCenterNodesFirstFloor.append(NnumMiddleStart)
	ops.geomTransf(coordTransf, pc*1000+NnumMiddleStart+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node2, NnumMiddleStart, 15, pc*1000+NnumMiddleStart+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],pc*1000+NnumMiddleStart+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node2, NnumMiddleStart], pc*1000+NnumMiddleStart+100000, 2)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],3,15,pc*1000+NnumMiddleStart+100000,'-iter', 1)	   

	NnumMiddleStart+=1	
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node1, (NnumMiddleStart-1), 15, pc*1000+(NnumMiddleStart-1)+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],pc*1000+(NnumMiddleStart-1)+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node1, (NnumMiddleStart-1)], pc*1000+(NnumMiddleStart-1)+100000, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],3,15,pc*1000+(NnumMiddleStart-1)+100000,'-iter', 1)	   
	NnumMiddleStart+=1   	
	ct+=1		  
	#vfo.plot_model()


BeamCenterNodesSecondFloor=[]

ct=0
NnumMiddleStart=20000000
for nodeCurrent in secondStoryColumnTopNodes[0:4]:
	node2=secondStoryColumnTopNodes[ct+1]
	node1=nodeCurrent
	print(nodeCurrent)
	node1Loc=ops.nodeCoord(node1)
	node2Loc=ops.nodeCoord(node2)
	middleLoc=(np.array(node2Loc)+np.array(node1Loc))/2

	beamNormal=np.cross(np.array(node2Loc)-np.array(node1Loc),[0,0,1])
    
	print(NnumMiddleStart,middleLoc,beamNormal) 
	ops.node(NnumMiddleStart, middleLoc[0],middleLoc[1],middleLoc[2])
	ops.mass(NnumMiddleStart,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])
	BeamCenterNodesSecondFloor.append(NnumMiddleStart)
	ops.geomTransf(coordTransf, pc*1000+NnumMiddleStart+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node2, NnumMiddleStart, 15, pc*1000+NnumMiddleStart+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],pc*1000+NnumMiddleStart+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node2, NnumMiddleStart], pc*1000+NnumMiddleStart+100000, 2)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],3,15,pc*1000+NnumMiddleStart+100000,'-iter', 1)	   



	NnumMiddleStart+=1	
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node1, (NnumMiddleStart-1), 15, pc*1000+(NnumMiddleStart-1)+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],pc*1000+(NnumMiddleStart-1)+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node1, (NnumMiddleStart-1)], pc*1000+(NnumMiddleStart-1)+100000, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],3,15,pc*1000+(NnumMiddleStart-1)+100000,'-iter', 1)	   
	NnumMiddleStart+=1   			  
	#vfo.plot_model()

    
# 	beamn=np.array(node2Loc)-np.array(node1Loc)
# 	GussLoc1=node1Loc+beamn*BeamConnLength+np.array([0,0,1])*ColumnConnectionLength
# 	GussLoc2=node1Loc+np.array([0,0,1])*ColumnConnectionLength
# 	GussLoc3=node1Loc
# 	GussLoc4=node1Loc+beamn*BeamConnLength
# 	CN1=60000000+NnumMiddleStart+Gct*4
# 	CN2=60000001+NnumMiddleStart+Gct*4
# 	CN3=60000002+NnumMiddleStart+Gct*4
# 	CN4=60000003+NnumMiddleStart+Gct*4
# 	ops.node(CN1, GussLoc1[0],GussLoc1[1],GussLoc1[2])
# 	ops.node(CN2, GussLoc2[0],GussLoc2[1],GussLoc2[2])
# 	ops.node(CN3, GussLoc3[0],GussLoc3[1],GussLoc3[2])
# 	ops.node(CN4, GussLoc4[0],GussLoc4[1],GussLoc4[2])
# 	eleTag=60000000+NnumMiddleStart+Gct
# 	eleNodes=[CN1,CN2,CN3,CN4]
# 	
# 	ops.element('ShellMITC4', eleTag, *eleNodes, secTag)
# 	print(NnumMiddleStart,middleLoc,beamNormal) 
    
# 	GL2CN=thirdStoryColumnBaseNodes[ct]
# 	ops.equalDOF(CN2, GL2CN, *dofs)
# 	Gct+=1
# 	beamn=np.array(node2Loc)-np.array(node1Loc)
# 	GussLoc1=node2Loc-beamn*BeamConnLength+np.array([0,0,1])*ColumnConnectionLength
# 	GussLoc2=node2Loc+np.array([0,0,1])*ColumnConnectionLength
# 	GussLoc3=node2Loc
# 	GussLoc4=node2Loc-beamn*BeamConnLength
# 	CN1=60000000+NnumMiddleStart+Gct*4
# 	CN2=60000001+NnumMiddleStart+Gct*4
# 	CN3=60000002+NnumMiddleStart+Gct*4
# 	CN4=60000003+NnumMiddleStart+Gct*4
# 	ops.node(CN1, GussLoc1[0],GussLoc1[1],GussLoc1[2])
# 	ops.node(CN2, GussLoc2[0],GussLoc2[1],GussLoc2[2])
# 	ops.node(CN3, GussLoc3[0],GussLoc3[1],GussLoc3[2])
# 	ops.node(CN4, GussLoc4[0],GussLoc4[1],GussLoc4[2])
# 	eleTag=60000000+NnumMiddleStart+Gct
# 	eleNodes=[CN1,CN2,CN3,CN4]
    
# 	ops.element('ShellMITC4', eleTag, *eleNodes, secTag)
# 	print(NnumMiddleStart,middleLoc,beamNormal) 
    
# 	GL2CN=thirdStoryColumnBaseNodes[ct+1]
# 	ops.equalDOF(CN2, GL2CN, *dofs)
# 	Gct+=1
	ct+=1
    
BeamCenterNodesThirdFloor=[]

ct=0
NnumMiddleStart=30000000
for nodeCurrent in thirdStoryColumnTopNodes[0:4]:
	node2=thirdStoryColumnTopNodes[ct+1]
	node1=nodeCurrent
	print(nodeCurrent)
    
    
	node1Loc=ops.nodeCoord(node1)
	node2Loc=ops.nodeCoord(node2)
    
    
	middleLoc=(np.array(node2Loc)+np.array(node1Loc))/2

	beamNormal=np.cross(np.array(node2Loc)-np.array(node1Loc),[0,0,1])
    
	print(NnumMiddleStart,middleLoc,beamNormal) 
	ops.node(NnumMiddleStart, middleLoc[0],middleLoc[1],middleLoc[2])
	ops.mass(NnumMiddleStart,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])
	BeamCenterNodesThirdFloor.append(NnumMiddleStart)
	ops.geomTransf(coordTransf, pc*1000+NnumMiddleStart+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node2, NnumMiddleStart, 15, pc*1000+NnumMiddleStart+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],pc*1000+NnumMiddleStart+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node2, NnumMiddleStart], pc*1000+NnumMiddleStart+100000, 2)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],3,15,pc*1000+NnumMiddleStart+100000,'-iter', 1)	   

	NnumMiddleStart+=1	
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node1, (NnumMiddleStart-1), 15, pc*1000+(NnumMiddleStart-1)+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],pc*1000+(NnumMiddleStart-1)+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node1, (NnumMiddleStart-1)], pc*1000+(NnumMiddleStart-1)+100000, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],3,15,pc*1000+(NnumMiddleStart-1)+100000,'-iter', 1)	   
	NnumMiddleStart+=1   			  

	ct+=1 

BeamCenterNodesTopFloor=[]
ct=0
NnumMiddleStart=40000000
for nodeCurrent in topStoryColumnTopNodes[0:4]:
	node2=topStoryColumnTopNodes[ct+1]
	node1=nodeCurrent
	print(nodeCurrent)
	node1Loc=ops.nodeCoord(node1)
	node2Loc=ops.nodeCoord(node2)
	middleLoc=(np.array(node2Loc)+np.array(node1Loc))/2

	beamNormal=np.cross(np.array(node2Loc)-np.array(node1Loc),[0,0,1])
    
	print(NnumMiddleStart,middleLoc,beamNormal) 
	ops.node(NnumMiddleStart, middleLoc[0],middleLoc[1],middleLoc[2])
	ops.mass(NnumMiddleStart,*[nodalMass,nodalMass,nodalMass,nodRotMass,nodRotMass,nodRotMass])
	BeamCenterNodesTopFloor.append(NnumMiddleStart)
	ops.geomTransf(coordTransf, pc*1000+NnumMiddleStart+100000, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node2, NnumMiddleStart, 15, pc*1000+NnumMiddleStart+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],pc*1000+NnumMiddleStart+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node2, NnumMiddleStart], pc*1000+NnumMiddleStart+100000, 2)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node2, NnumMiddleStart],3,15,pc*1000+NnumMiddleStart+100000,'-iter', 1)	   



	NnumMiddleStart+=1	
	if beamType=='elastic':
		ops.element('elasticBeamColumn', NnumMiddleStart, node1, (NnumMiddleStart-1), 15, pc*1000+(NnumMiddleStart-1)+100000)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],pc*1000+(NnumMiddleStart-1)+100000, 2)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+NnumMiddleStart, *[node1, (NnumMiddleStart-1)], pc*1000+(NnumMiddleStart-1)+100000, 1)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', NnumMiddleStart, *[node1,(NnumMiddleStart-1)],3,15,pc*1000+(NnumMiddleStart-1)+100000,'-iter', 1)	   
	NnumMiddleStart+=1   			  
	#vfo.plot_model()
	ct+=1



print('I am here columns are made')
#vfo.plot_model()


Zbase=zBase
Ztop=z1stFloor


pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]

if useSlabs==1:
	# Breakaway first STORY SLAB
	##################################################################################################################################################################
	zBumpDebug=0.0
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


	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]


	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=1
	startEle=1

	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)
	firstSlabCornerNodeNums=[startNode,startNode+numX, -1 + startNode+(numX+1)*(numY+1),-1+startNode+(numX+1)*(numY+1) - numX]
	firstSlabEdgeMidNodeNums=[startNode+numY/2,startNode+(((numX+1)*(numY+1))-1)-(numX/2*(numY+1)), startNode+(((numX+1)*(numY+1))-1)-numY/2,startNode+(((numX+1)*(numY+1))-1)-((numX/2+1)*(numY+1))+1]

	SlabNodes=np.linspace(startNode,startNode+((numX+1)*(numY+1))-1,((numX+1)*(numY+1)))
	for nodeNum in SlabNodes:
		FOAMySeesInstance.coupledNodes.append(int(nodeNum))



#vfo.plot_model()
Zbase=z1stFloor
#Ztop=1.8288
Ztop=z2ndFloor


if useSlabs==1:
	# 2nd STORY SLAB
	##################################################################################################################################################################

	zBumpDebug=0.0
	Ztop+=zBumpDebug
	pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]



	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]


	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=1000
	startEle=1000

	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)
	secondSlabCornerNodeNums=[startNode,startNode+numX, -1 + startNode+(numX+1)*(numY+1),-1+startNode+(numX+1)*(numY+1) - numX]
	secondSlabEdgeMidNodeNums=[startNode+numY/2,startNode+(((numX+1)*(numY+1))-1)-(numX/2*(numY+1)), startNode+(((numX+1)*(numY+1))-1)-numY/2,startNode+(((numX+1)*(numY+1))-1)-((numX/2+1)*(numY+1))+1]



#vfo.plot_model()	
Zbase=z2ndFloor
#Ztop=1.8288
Ztop=z3rdFloor


if useSlabs==1:
	# Third STORY SLAB
	##################################################################################################################################################################

	zBumpDebug=0.0
	Ztop+=zBumpDebug

	pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]




	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]



	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=222000
	startEle=222000

	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)
	thirdSlabCornerNodeNums=[startNode,startNode+numX, -1 + startNode+(numX+1)*(numY+1),-1+startNode+(numX+1)*(numY+1) - numX]
	thirdSlabEdgeMidNodeNums=[startNode+numY/2,startNode+(((numX+1)*(numY+1))-1)-(numX/2*(numY+1)), startNode+(((numX+1)*(numY+1))-1)-numY/2,startNode+(((numX+1)*(numY+1))-1)-((numX/2+1)*(numY+1))+1]

ct=0

    
Zbase=z3rdFloor
#Ztop=1.8288
Ztop=z4thFloor


if useSlabs==1:
	# Top STORY SLAB
	##################################################################################################################################################################

	zBumpDebug=0.0
	Ztop+=zBumpDebug

	pileLocs=[[Xmin,Ymin,Ztop],[Xmin,Ymax,Ztop],[Xmax,Ymax,Ztop],[Xmax,Ymin,Ztop]]




	coooords=[1,pileLocs[0][0],pileLocs[0][1],pileLocs[0][2],2,pileLocs[1][0],pileLocs[1][1],pileLocs[1][2],3,pileLocs[2][0],pileLocs[2][1],pileLocs[2][2],4, pileLocs[3][0],pileLocs[3][1],pileLocs[3][2]]



	# block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *crds)
	eleType='shell'
	eleArgs=[100]
	startNode=1222000
	startEle=1222000
	ops.block2D(numX, numY, startNode, startEle, eleType, *eleArgs, *coooords)
	topSlabCornerNodeNums=[startNode,startNode+numX, -1 + startNode+(numX+1)*(numY+1),-1+startNode+(numX+1)*(numY+1) - numX]
	topSlabEdgeMidNodeNums=[startNode+numY/2,startNode+(((numX+1)*(numY+1))-1)-(numX/2*(numY+1)), startNode+(((numX+1)*(numY+1))-1)-numY/2,startNode+(((numX+1)*(numY+1))-1)-((numX/2+1)*(numY+1))+1]


ct=0
for corner in firstSlabCornerNodeNums:
	columnTopNode=firstStoryColumnNodes[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1

ct=0
for corner in secondSlabCornerNodeNums:
	columnTopNode=secondStoryColumnTopNodes[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1
    
ct=0
for corner in thirdSlabCornerNodeNums:
	columnTopNode=thirdStoryColumnTopNodes[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1

ct=0
for corner in topSlabCornerNodeNums:
	columnTopNode=topStoryColumnTopNodes[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1


ct=0
for corner in firstSlabEdgeMidNodeNums:
	columnTopNode=BeamCenterNodesFirstFloor[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1

ct=0
for corner in secondSlabEdgeMidNodeNums:
	columnTopNode=BeamCenterNodesSecondFloor[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1
    
ct=0
for corner in thirdSlabEdgeMidNodeNums:
	columnTopNode=BeamCenterNodesThirdFloor[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1

ct=0
for corner in topSlabEdgeMidNodeNums:
	columnTopNode=BeamCenterNodesTopFloor[ct]
	beamStartNode=corner
	rNodeTag=columnTopNode
	cNodeTag=beamStartNode
	if constraintType=='EQDOF': 
		ops.equalDOF(rNodeTag, cNodeTag, *dofs)
	elif constraintType=='RIGIDLINK':	
		ops.rigidLink('bar', rNodeTag, cNodeTag)
	elif constraintType=='EQDOFMIXED':
		ops.equalDOF_Mixed(rNodeTag, cNodeTag, numDOF, *rcdofs)
	ct+=1



#vfo.plot_model()
# if useBeams==1:

    
if noBraces==1:
	pass
else:
	############### Chevron braces #####################


	# corner nodes nElemBOTTOM  x 1000+nElemBOTTOM  y  2000+nElemBOTTOM  x  3000+nElemBOTTOM   y   nElemBOTTOM

	# mid nodes			   20000+numX/2	 21000+numX/2	 22000+numX/2		  23000+numX/2
			
	pc+=1

	CT=101010101
	eleNo=CT
	beamNormal=[1,0,0]
	nodeNum+=1	
	ops.beamIntegration('Lobatto', 16, secTag, 2)
	Nodes=[25000,BeamCenterNodesThirdFloor[0]]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
		
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1
    
	nodeNum+=1
	CT=101010102
	eleNo=CT
	beamNormal=[1,0,0]
	Nodes=[BeamCenterNodesThirdFloor[0],26000]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
		
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1

	nodeNum+=1
	CT=101010105
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[26000,BeamCenterNodesThirdFloor[1]]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
		
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1
    
	nodeNum+=1
	CT=101010106
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[BeamCenterNodesThirdFloor[1],27000]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1
    


	nodeNum+=1
	CT=101010103
	eleNo=CT
	beamNormal=[1,0,0]
	Nodes=[27000,BeamCenterNodesThirdFloor[2]]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1
    

	nodeNum+=1
	CT=101010104
	eleNo=CT
	beamNormal=[1,0,0]
	Nodes=[BeamCenterNodesThirdFloor[2],28000]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1
    
	nodeNum+=1
	CT=101010107
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[28000,BeamCenterNodesThirdFloor[3]]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1

	nodeNum+=1
	CT=101010108
	eleNo=CT
	beamNormal=[0,1,0]
	Nodes=[BeamCenterNodesThirdFloor[3],25000]
	TransfTag=pc*1000+nodeNum+100000
	ops.geomTransf(coordTransf, TransfTag, beamNormal[0],beamNormal[1],beamNormal[2])
	if beamType=='elastic':
		ops.element('elasticBeamColumn', nodeNum, Nodes[0], Nodes[1], 101, TransfTag)  		
	elif beamType=='dispBeamCol':
		ops.element('dispBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		 #, '-cMass', '-mass', mass=0.0)
	elif beamType=='forceBeamCol':
		ops.element('forceBeamColumn', pc*1000+nodeNum, *[Nodes[0], Nodes[1]],TransfTag, 3)		
	elif beamType=='NLBeamCol':
		ops.element('nonlinearBeamColumn', nodeNum, *[Nodes[0], Nodes[1]],3,101,TransfTag,'-iter', 1)	   
		pc+=1
    


#vfo.plot_model()	
print(FOAMySeesInstance.coupledNodes)

if BaseIsolated==1:
	###############################################################################################################################################
	#   BASE ISOLATOR ELEMENTS #
	###############################################################################################################################################



	ops.node(9000000, *[Xmin,Ymin,zBase])
	ops.node(9001000, *[Xmin,Ymax,zBase])
	ops.node(9002000, *[Xmax,Ymax,zBase])
	ops.node(9003000, *[Xmax,Ymin,zBase])


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
ops.fixZ(1.0,*[1,1,1,1,1,1])

ops.recorder('Node', '-file', 'frontRightCFTReactionAtFlumeFloor.out','-time', '-node', 1000000, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'frontLeftCFTReactionAtFlumeFloor.out','-time', '-node',1001000, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'backLeftCFTReactionAtFlumeFloor.out','-time', '-node', 1002000, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'backRightCFTReactionAtFlumeFloor.out','-time', '-node',1003000, '-dof', 1,2,3,4,5,6, 'reaction')


ops.recorder('Node', '-file', 'frontRightCFTReaction.out','-time', '-node', 123456780, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'frontLeftCFTReaction.out','-time', '-node',123456781, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'backLeftCFTReaction.out','-time', '-node', 123456782, '-dof', 1,2,3,4,5,6, 'reaction')
ops.recorder('Node', '-file', 'backRightCFTReaction.out','-time', '-node',123456783, '-dof', 1,2,3,4,5,6, 'reaction')
res=['disp','vel','accel','incrDisp','reaction','pressure','unbalancedLoad','mass']

ops.recorder('PVD', 'output', '-precision', 4, '-dT', .05, *res)
beep=0
for SGE in StrainGaugeElements:
    beep+=1
    ops.recorder('Element','-ele',SGE,'-file',str(beep)+'fiber.out','section','fiber',OuterRad,0,600,'stressStrain')

for soiln in soilnodes:
    
    ops.recorder('Node', '-file', 'SoilNode'+str(soiln)+'Reaction.out','-time', '-node',soiln, '-dof', 1,2,3,4,5,6, 'reaction')



# GMfile='./GroundMotion.acc'
# # Uniform EXCITATION: acceleration input
# IDloadTag = 400            # load tag
# dt = 0.001            # time step for input ground motion
# GMfact = 1           # data in input file is in g Unifts -- ACCELERATION TH
maxNumIter = 10
# GMdirection=1
Tol=1e-3

# TSTAG=2

# ops.timeSeries('Path', TSTAG, '-dt', dt, '-filePath', GMfile, '-factor', GMfact, '-useLast')

# ops.pattern('UniformExcitation', IDloadTag, GMdirection, '-accel', 2) 


# ops.constraints('Transformation')
# ops.numberer('Plain')
# ops.system('BandGeneral')
ops.test('EnergyIncr', Tol, maxNumIter)
# ops.algorithm('ModifiedNewton')
# NewmarkGamma = 0.5
# NewmarkBeta = 0.25
# ops.integrator('Newmark', NewmarkGamma, NewmarkBeta)
# ops.analysis('VariableTransient')
# DtAnalysis = 0.001
# TmaxAnalysis = 10
# Nsteps =  int(TmaxAnalysis/ DtAnalysis)
# ok=1
# # for i in test:
# ops.algorithm('KrylovNewton')


# ok = ops.analyze(Nsteps, DtAnalysis,DtAnalysis/1000,DtAnalysis,100)    
# ops.remove('loadPattern',IDloadTag)
# ops.remove('timeSeries', 2) 
	
	