from odbAccess import *
from abaqusConstants import *
from numpy import *
import numpy as np
import math

#odb data
#pathodb='D:\DBGuan\Growth_Github\Growth-ConstrainMatrix\Growth_Update\Eccentric_Plastic_Original/'
#pathfile=pathodb
#
file='RBM_updata.odb'
odb=openOdb(file)
assembly=odb.rootAssembly
# instance name
theinst=assembly.instances['PART-1_1']
# node set name
nodest=theinst.nodeSets['BOCAP']
# element set name
elementst=theinst.elementSets['EALL']
#step name
step1=odb.steps['Relax']
theframe=step1.frames

#frame number; -1 is index of last frame
#extract out displacement values
TempField=theframe[-1].fieldOutputs['U'] 
ns2disp=TempField.getSubset(region=nodest)
ns2value=ns2disp.values

#create write out data file
#output=open(pathfile+'node.txt','w')
#index=-1
#for s in ns2value:
#    index=index+1
#    dispcomp=s.data
#    ndcoord=s.instance.nodes[index].coordinates
#    ndx=ndcoord[0]+dispcomp[0]
#    ndy=ndcoord[1]+dispcomp[1]
#    ndz=ndcoord[2]+dispcomp[2]
#    ndID=s.instance.nodes[index].label
#    output.write('%i,\t %14.10f,\t %14.10f,\t %14.10f\n' %(ndID, ndx,ndy,ndz))

#output.close()

output8=open('Ndata.txt','w')
#output8.write('%i\n' %(len(ns2value)))
index=-1
for s in ns2value:
    index=index+1
    dispcomp=s.data
    ndcoord=s.instance.nodes[index].coordinates
    ndx=ndcoord[0]
    ndy=ndcoord[1]
    ndz=ndcoord[2]
    ndID=s.instance.nodes[index].label
    output8.write('%i,\t %14.10f,\t %14.10f,\t %14.10f\n' %(ndID, ndx,ndy,ndz))

output8.close()


output7=open('Udata.txt','w')
#output7.write('%i\n' %(len(ns2value)))
index=-1
for s in ns2value:
    index=index+1
    dispcomp=s.data
    ndcoord=s.instance.nodes[index].coordinates
    ndx=dispcomp[0]
    ndy=dispcomp[1]
    ndz=dispcomp[2]
    ndID=s.instance.nodes[index].label
    output7.write('%i,\t %14.10f,\t %14.10f,\t %14.10f\n' %(ndID, ndx,ndy,ndz))

output7.close()


step3=odb.steps['Recovery']
theframe3=step3.frames

SDVField1=theframe3[-2].fieldOutputs['SDV1'] 
#extract out displacement values
sdv1=SDVField1.getSubset(region=elementst)
sdv2value1=sdv1.values

SDVField2=theframe3[-2].fieldOutputs['SDV2'] 
#extract out displacement values
sdv2=SDVField2.getSubset(region=elementst)
sdv2value2=sdv2.values

SDVField3=theframe3[-2].fieldOutputs['SDV3'] 
#extract out displacement values
sdv3=SDVField3.getSubset(region=elementst)
sdv2value3=sdv3.values

#create write out data file
output3=open('Stretch.txt','w')
#output3.write('%i\n' %(len(sdv2value)))
index=-1
for s in sdv2value1:
    index=index+1
    sdvalue1=s.data
    ndID=s.elementLabel
    sdvalue2=sdv2value2[index].data	
    sdvalue3=sdv2value3[index].data
    output3.write('%i,\t %14.10f,\t %14.10f,\t %14.10f\n' %(ndID, sdvalue1, sdvalue2, sdvalue3))

output3.close()

#SField=theframe[12].fieldOutputs['S']
#sdv5=SField.getSubset(region=elementst)
#sdv5value=sdv5.values
#output5=open(pathfile+'Stress.txt','a')
#index=-1
#stre1=0.0
#stre2=0.0
#stre3=0.0
#stre4=0.0
#stre5=0.0
#stre6=0.0
#for s in sdv5value:
#    index=index+1
#    sdvalue=s.data
#    ndID=s.elementLabel
#    stre1=stre1+sdvalue[0]
#    stre2=stre2+sdvalue[1]
#    stre3=stre3+sdvalue[2]
#    stre4=stre4+sdvalue[3]
#    stre5=stre5+sdvalue[4]
#    stre6=stre6+sdvalue[5]

	
#output5.write('%14.10f,\t %14.10f,\t %14.10f,\t %14.10f,\t %14.10f,\t %14.10f\n' \
#    %(stre1/len(sdv5value),stre2/len(sdv5value),stre3/len(sdv5value), \
#      stre4/len(sdv5value),stre5/len(sdv5value),stre6/len(sdv5value)))

###################################################
# ACTIVE SYSTOLE
step2=odb.steps['Beat']
theframe2=step2.frames

SField=theframe2[-1].fieldOutputs['SDV4'] 
sdv5=SField.getSubset(region=elementst)
sdv5value=sdv5.values
output5=open('Stress.txt','w')
index=-1
for s in sdv5value:
    index=index+1
    sdvalue=s.data
    ndID=s.elementLabel
    output5.write('%i,\t %14.10f\n' %(ndID, sdvalue))

output5.close()

v=step2.historyRegions
#sd2=v['Element PART-1_1.21925 Int Point 1'].historyOutputs['SDV2']
#output2=open(pathfile+'SDV2.txt','w')
#for time, value in sd2.data:
#    output2.write('%14.10f, \t%14.10f\n' %(time, value))
#output2.close()

#sd1=v['Element VENTRICLES-1.3012 Int Point 1'].historyOutputs['SDV1']
#output1=open(pathfile+'SDV1.txt','a')
#for time, value in sd1.data:
#    output1.write('%14.10f, \t%14.10f\n' %(time, value))
#output1.close()

sd4=v['Node ASSEMBLY.1'].historyOutputs['CVOL']
output4=open('CVOLV.txt','a')
for time, value in sd4.data:
    output4.write('%14.10f, \t%14.10f\n' %(time, value))
output4.close()


sd4=v['Node ASSEMBLY.1'].historyOutputs['PCAV']
output41=open('PCALV.txt','a')
for time, value in sd4.data:
    output41.write('%14.10f, \t%14.10f\n' %(time, value))
output41.close()
#write out the sdv values for every element

v2=step3.historyRegions

sd4=v2['Node ASSEMBLY.1'].historyOutputs['CVOL']
output4=open('CVOLV.txt','a')
for time, value in sd4.data:
    output4.write('%14.10f, \t%14.10f\n' %(time, value))
output4.close()


sd4=v2['Node ASSEMBLY.1'].historyOutputs['PCAV']
output41=open('PCALV.txt','a')
for time, value in sd4.data:
    output41.write('%14.10f, \t%14.10f\n' %(time, value))
output41.close()
#write out the sdv values for every element


odb.close()



