from ij import WindowManager, IJ,Prefs,ImagePlus
from ij.measure import ResultsTable, Measurements
from ij.plugin import ZProjector,RGBStackMerge, ImageCalculator
from ij.plugin.frame import RoiManager
from java.awt import Color
from ij.process.LUT import createLutFromColor
from ij.gui import WaitForUserDialog,Plot

import math,os,re,sys,time
from os.path import exists

# function to get orientation of the image manually


def getAngle(imp=None,nRep=3):
	if (imp==None):
		imp=IJ.getImage()
	angs=[]
	IJ.setTool("line")
	for i in range(nRep):
		msgN="Select the apical-basal axis with Line tool! ("+str(i+1)+"/"+str(nRep)+")"	 
		msg=WaitForUserDialog(msgN)
		msg.show()
		roi=imp.getRoi()
		# print(roi==None)
		
		while (WindowManager.getWindowCount()>0) and (roi==None or roi.isLine()==False):
			msg=WaitForUserDialog("ROI type wrong! Select the axis with Line tool!")
			msg.show()
			roi=imp.getRoi()
		if (not roi==None):
			angs.append(roi.getAngle())
			imp.deleteRoi()
	if (len(angs)>0):
		ang_ave=sum(angs)/len(angs)
		print("the angle is:"+str(ang_ave)+" from "+str(len(angs))+" replicates")
		rt_main.addValue("Angle", ang_ave)
		rt_main.addValue("nAng", len(angs))
		return(ang_ave)
		
	else:
		return 0

# function to remove a folder recursively

def delFiles(dir): 
	print("Deleting directory: "+dir)
	if dir.endswith("/") == False: dir = dir+"/"
	list = os.listdir(dir)
	# print("number of item in list: % s" % len(list))
	if (len(list)==0):
		# print("Deleting empty directory: "+dir)
		os.rmdir(dir);
	else:	
		for f in list: 
			print("now comes to "+f)
			if (os.path.isfile(dir+"/"+f)):
				print("Deleting file: "+dir+"/"+f)
				os.remove(dir+"/"+f)
			else:
				# print("Deleting file:"+dir+"/"+f)
				delFiles(dir +"/"+f)
	print("% s removed successfully in function" % dir) 
	return 0		

#Function to filter out certain color from a RGB image
def ColorFilter(imp=None,col=[255,0,0]):
	if imp==None:
		imp=IJ.getImage()
	imp_title=imp.getTitle()
	
	imp1=imp.duplicate()
	imp1.setTitle("dup")
	IJ.run(imp1, "RGB Stack", "");
	for i in range(len(col)):
		minT=max(0,col[i]-3)
		maxT=min(255,col[i]+3)
		print("managing channel "+str(i))
		imp1.setC(i+1)
		IJ.setRawThreshold(imp1, minT, maxT,"")
		IJ.run(imp1, "Convert to Mask", "method=Default background=Light only")
	imp2=ZProjector.run(imp1,"min")
	imp2.setTitle(imp_title)
	return imp2


def GetMask(imp=None):
	if (imp==None):
		imp=IJ.getImage()
	imp1=imp.duplicate()
	IJ.run(imp1, "Gaussian Blur...", "sigma=5")
	imp1.getProcessor().setThreshold(0, 10, imp.getProcessor().NO_LUT_UPDATE)
	IJ.run(imp1, "Convert to Mask", "")
	IJ.run(imp1, "Invert", "")
	IJ.run(imp1, "Options...", "iterations=5 count=1 black do=Erode")
	IJ.run(imp1, "Analyze Particles...", "size=1600-Infinity show=Masks clear in_situ");
	IJ.run(imp1, "Invert", "")
	# imp1.show()
	return(imp1)

def GetOutline(imp=None):
	if (imp==None):
		imp=IJ.getImage()
	
	# imp.show()
	imp2=imp.duplicate()
	IJ.run(imp2, "Convert to Mask", "")
	#imp2=GetMask(imp1)
	imp2.setRoi(0,0,10,10)
	mes=imp2.getStatistics().mean
	imp2.deleteRoi()
	if (mes>200): #white background
		IJ.run(imp2, "Invert", "")
	
	# imp2.show()
	# IJ.run(imp2, "Select None", "")
	IJ.run(imp2, "Outline", "")
	mes=imp2.getStatistics().mean
	if (mes>200): #white background
		IJ.run(imp2, "Invert", "")
	IJ.run(imp2,"Grays","")
	return(imp2)
	
	
def colorPanel():
	imp = IJ.createImage("Test1", "RGB white", 1280, 1280, 1)
	for i in range(256):
		for j in range(256):
			r=i/float(255)
			g=j/float(255)
			b=1.5-r-g
			if b<=1 and b>=0:
				# print((r,g,b))
				cur_col=Color(r,g,b)
				imp.setColor(cur_col)
				imp.getProcessor().fillRect(i*5, j*5, 5, 5)
	
	print("Color panel generated!")
	imp.show()
	return(0)
	
# colorPanel()
# sys.exit()		
blue=[ 0,153,255 ]# Blue
yellow=[ 255,153,0 ] #Yellow
pinkred=[ 255,40,89 ] #pinkred
purple=[ 119,12,252 ] #purple
green=[ 145,238,0 ] #green
turqoise=[ 9,193,181 ] #turqoise
orange=[ 246,132,5 ] #orange
	
def ColorScale(cS=blue,cE=orange,colLen=11):
	
	Rs=[]
	Gs=[]
	Bs=[]
	
	for i in range(colLen):
		# print(rgb)
		incre=i/float(colLen-1)/float(255)
		# print(type(incre))
		Rs.append(incre*(cE[0]-cS[0])+cS[0]/float(255))
		Gs.append(incre*(cE[1]-cS[1])+cS[1]/float(255))
		Bs.append(incre*(cE[2]-cS[2])+cS[2]/float(255))
	return [Rs,Gs,Bs]

def DrawScale(CC=ColorScale(orange,blue,11),w=400,h=100):
	# print("heiehiehi")
	if (not len(CC)==3):
		print("Dimension of the color code array wrong!")
		return 0
	nCol=len(CC[0]) # number of total colors
	
	#Some dimensions
	hBor=10
	wBor=25
	hLine=20
	hLetter=15
	
	imp = IJ.createImage("ColorScale1", "RGB white", w+wBor*2, h+hLine+hBor*2+hLetter*2, 1);
	wblock=w/nCol
	hblock=h
	
	rt=ResultsTable()
	
	# Draw color blocks
	for c in range(nCol):
		cur_col=Color(CC[0][c],CC[1][c],CC[2][c])
		imp.setColor(cur_col)
		imp.getProcessor().fillRect(c*wblock+25, 10, wblock, hblock)
		rt.incrementCounter()
		rt.addValue("Red",CC[0][c])
		rt.addValue("Green",CC[1][c])
		rt.addValue("Blue",CC[2][c])
	
	#Draw values
	posL=25+int(0.5*wblock)
	posM=25+int(0.5*wblock*nCol)
	posR=25+int(wblock*nCol-0.5*wblock)
	
	imp.setColor(Color.black)
	imp.getProcessor().drawLine(posL,10+hblock,posL,10+hLine+hblock)
	#imp.getProcessor().drawLine(posM,10+hblock,posM,10+hLine+hblock)
	imp.getProcessor().drawLine(posR,10+hblock,posR,10+hLine+hblock)
	
	#Add labels
	imp.getProcessor().setJustification(1)
	imp.getProcessor().drawString("0",posL,10+hLine+hblock+hLetter)
	#imp.getProcessor().drawString(,posM,10+hLine+hblock+hLetter)
	imp.getProcessor().drawString(str(int((nCol-1)*5)),posR,10+hLine+hblock+hLetter)
	imp.getProcessor().drawString("Distance to the dorsal surface (micron)",posM,10+hLine+hblock+hLetter*2)
		
	# imp.show()
	# rt.show("Color Palette")
	return imp
	

def median(lst,per):
    n = len(lst)
    s = sorted(lst)
    if n>1:
    	hf=min(n-1,int(round(n*per)))
    	return s[hf]
    else:
    	return 0
    	

def mean(lst):
    n = len(lst)
    s = sum(lst)
    return float(s)/n if n else None
    
def std(lst):
	n=len(lst)
	if n:
		m=mean(lst)
		sumsq = sum((v-m) ** 2 for v in lst)
		if (n>1):
	   		return (sumsq / (n-1)) ** 0.5
	   	else:
	   		return 0
	else:
		return None

#
#def FindMid(imp=None,outdir=IJ.getDirectory("file")):
#	if imp==None: imp=IJ.getImage()
#	imp_title=imp.getTitle()
#	imp_title=re.sub(".*[\\\/]","",imp_title)
#	
#	imp_bi=imp.duplicate()
#	#Reduce the size of image for fast looping
#	imp_bi=imp_bi.resize(int(0.5*imp_bi.width), int(0.5*imp_bi.height), 1, "bilinear")
#	
#	pd=imp_bi.getCalibration().pixelDepth
#	
#	#Reslice to avoid insufficient stacking
#	rsd=1
#	IJ.run(imp_bi, "Reslice Z", "new="+str(rsd))
#	imp_bi=IJ.getImage()
#	imp_bi.setTitle("imp_bi")
#	
#	#Make binary
#	IJ.run(imp_bi, "Make Binary", "method=IsoData background=Dark calculate black")
#	#imp_bi.show()
#	
#	nS=imp_bi.getNSlices()
#	ht=imp_bi.height
#	wd=imp_bi.width
#	ip=imp_bi.getProcessor()
#	
#	MESs=[0]*nS
#	Xsd=[0]*nS
#	#Calculate standard diviation of X positions
#	for fi in range(nS):
#		imp_bi.setSliceWithoutUpdate(fi+1)
#		stat=ip.getStatistics()
#		MESs[fi]=stat.mean
##		Xs=[]
##		for fx in range(wd):
##			for fy in range(ht):
##				value=ip.getPixel(fx,fy)
##				if value>200:
##					Xs.append(fx)
##		if (len(Xs)>0):
##			mx=mean(Xs)
##			Xdist=[abs(x-mx) for x in Xs]
##			Xsd[fi]=median(Xdist,0.75)
#		
#	
#	#Averaging neighboring values
#	mR=3 #int(5/rsd)
#	avgMES=[0]*nS
#	#avgXsd=[]
#	for fi in range(nS):
#		left=max(0, fi-mR)
#		right=min(fi+mR,nS-1)+1
#		mMES=mean(MESs[left:right])
##		mXsd=mean(Xsd[left:right])
#		#avgMES[fi]=mMES*mXsd
#		avgMES[fi]=mMES#mXsd
#		#avgXsd.append()
#	
#	#Plot for Z-profile
#	plot1=Plot("Averge Z-plot","new Slice",str(2*mR+1)+" Slice average")
#	plot1.addPoints(range(1,nS+1),avgMES,Plot.CROSS)
#	pw=plot1.show()
#
#	#Find local minima in averages
#	minSlc=[]
#	maxSlc=[]
#	stat=imp_bi.getStatistics()
#	nblank=[k+1 for k in range(len(avgMES)) if avgMES[k]>0.01]
#	rL0=min(nblank)
#	rR0=max(nblank)
#	print(["nblank",rL0,rR0])
#	
#	for fi in range(rL0-1,rR0):
#		left=max(0, fi-mR)
#		right=min(fi+mR,nS-1)+1
#		if stat.mean <50 and avgMES[fi]==min(avgMES[left:right]): # a white image
#				minSlc.append(fi+1)
#		if stat.mean <50 and avgMES[fi]==max(avgMES[left:right]): # a white image
#				maxSlc.append(fi+1)
#	
#	print(["minSlc",minSlc])
#	print(["maxSlc",maxSlc])
#
#	plot1.setColor(Color(255,153,0))
#	plot1.addPoints(maxSlc,[avgMES[x-1] for x in maxSlc],Plot.CIRCLE)
#	pw=plot1.update()
#
#	#Expand the range considering fluctuation of values
#	#Considering the values between two peaks
#	if (len(maxSlc)>1 and (max(maxSlc)-min(maxSlc)>pd/rsd)): #at least two local maxima
#		rL=min(maxSlc)
#		rR=max(maxSlc)
#	else: #only one local maxima detected
#		if (min(maxSlc)-1)<(nS-max(maxSlc)): #peak close to left end
#			rL=min(maxSlc)
#			rR=nS
#		else:
#			rL=1
#			rR=max(maxSlc)
#	#Mark the range between two peaks
#	plot1.setColor(Color(0,0,255))
#	plot1.addPoints(range(rL,rR+1),[avgMES[x-1] for x in range(rL,rR+1)],Plot.CROSS)
#
#	#Range of fluctuation
#	maxMES=max(avgMES)
#	minSlc1=[]
#	fluct=0.025 # fluctuation of 5 pixels
#	for s in minSlc:
#		val=avgMES[s-1]
#		minV=val-fluct*val
#		maxV=val+fluct*val
#
#		d=0
#		while(s-d>rL-1):
#			if minV<avgMES[s-d-1]<maxV:
#				minSlc1.append(s-d)
#			d=d+1
#		d=1
#		while(s+d<rR-1):
#			if minV<avgMES[s+d-1]<maxV:
#				minSlc1.append(s+d)
#			d=d+1
#			
#	#Select the local min/max most close to the middle of stack
#	print(["new minSlc",minSlc1])
#	#Mark the range of Slc1
#	plot1.setColor(Color(0,255,0))
#	plot1.addPoints(minSlc1,[avgMES[x-1] for x in minSlc1],Plot.CROSS)
#
#	
#	midStack=int(round(float(rL0+rR0)*0.5))
#	print("midStack in resliced stack:"+str(midStack))
#	if (len(minSlc)>0):
#		DifSlc=[abs(x-midStack) for x in minSlc1]
#		midSlc=minSlc1[DifSlc.index(min(DifSlc))]
#	else:
#		midSlc=midStack
#		
#	print("Middle slice found in resliced stack at: "+str(midSlc))	
#	midSlc_ori=int(round((midSlc-1)*rsd/pd+1))
#	print("Middle slice found at: "+str(midSlc_ori))
#	
#	#Save the plot of average Z-profiles
#	#print(len(avgMES))
#	#print(avgMES)
#	plot1.setColor(Color(255,0,0))
#	plot1.addPoints([midStack],avgMES[(midStack-1):midStack],Plot.CIRCLE)
#	pw=plot1.update()
#	pl=IJ.getImage()
#	out=os.path.join(outdir,re.sub("\..*","-plot.jpg",imp_title))
#	print("Saving to: "+out)
#	IJ.saveAs(pl, "Jpeg", out)
#	#pw.dispose()
#
##	#Save the plot of Z-profiles
##	plot=Plot("Z-plot","Slice","Std(X)")
##	plot.addPoints(range(1,nS+1),MESs,Plot.CROSS)
##	plot.setColor(Color(255,0,0))
##	pw=plot.show()
##	pl=IJ.getImage()
##	out=os.path.join(outdir,re.sub("\..*","-plot.jpg",imp_title))
##	print("Saving to: "+out)
##	IJ.saveAs(pl, "Jpeg", out)
#	
#	return midSlc_ori
	
def ColorCode(imp=None,start=1,end=1,cS=blue,cE=orange,limit=20):
	if imp==None: imp=IJ.getImage()
	imp_title=imp.getTitle()
	
	lut_ori=imp.getProcessor().getLut()
	nS=imp.getNSlices()
	pd=imp.getCalibration().pixelDepth
	midSlc=0.5*(start+end)
	
		
	newStack = IJ.createImage("HyperStack", "RGB color-mode", imp.width, imp.height, 1, 1, nS)
	newStack.copyScale(imp)
	LB=int(max(1,start,math.ceil(((midSlc-1)*pd-limit)/pd+1)))
	RB=int(min(nS,end,math.floor(((midSlc-1)*pd+limit)/pd+1)))
	colLen=RB-LB+1
	CC=ColorScale(cS,cE,colLen)
	ds=DrawScale(CC,400,50)
	IJ.run(ds, "Multiply...", "value=1.5 stack")
	out = os.path.join(OutLineDir,basename+"-colors.tif")
	print("Saving to: "+out)
	IJ.saveAs(ds, "Tiff", out)
	
	for i in range(nS):
		#Convert to sequence in color code
		colorcode=i-LB
		if (i<=LB): colorcode=0
		if (i>=RB): colorcode=colLen-1
		
		lut=createLutFromColor(Color(CC[0][colorcode],CC[1][colorcode],CC[2][colorcode]))
		imp.setSlice(i+1)
		imp.setLut(lut)
		IJ.run(imp, "Select All", "")
		imp.copy();
		newStack.setPosition(1, 1, i + 1)
		newStack.paste()
	IJ.run(imp, "Select None", "")
	imp.setLut(lut_ori)
	IJ.run(newStack, "Select None", "")
	newStack.setTitle("CC_"+imp_title)
	# newStack.show()
	# newMerge = ZProjector.run(newStack,"max");			
	# newMerge.show()
	IJ.run(newStack, "Multiply...", "value=1.5 stack")
	return newStack	

# Test the function
# degree=getOrientation()
# print("degree="+str(degree))
# imp = WindowManager.getCurrentImage()
# IJ.run(imp, "Rotate... ", "angle="+str(degree-90)+" grid=1 interpolation=Bilinear enlarge");

#load files from folder
HomeDir = "Path/to/Original_stacks"

print(exists(HomeDir))

if exists(HomeDir)==False:
	print("Cannot find your files! Check HomeDir")
	sys.exit()

OutDir = HomeDir+"-Results"

# Remove old outputs
if exists(OutDir):
	print("Wants to delete "+OutDir);
	try: 
		delFiles(OutDir) 
		print("% s removed successfully" % OutDir) 
	except OSError as error: 
		print(error) 
		print("File path can not be removed") 
# Create new output folder
if (exists(OutDir)==0):
	try:
		os.mkdir(OutDir)
		print("% s created successfully" % OutDir) 
	except OSError as error: 
		print(error) 
		print("% s can not be created" % OutDir)
# Create new output folder for coordinates
CoOutDir=OutDir+"/Coordinates"
if (exists(CoOutDir)==0):
	try:
		os.mkdir(CoOutDir)
		print("% s created successfully" % OutDir) 
	except OSError as error: 
		print(error) 
		print("% s can not be created" % OutDir)
		
# Create new output folder for stardist outputs
StOutDir=OutDir+"/StarDist"
if (exists(StOutDir)==0):
	try:
		os.mkdir(StOutDir)
		print("% s created successfully" % OutDir) 
	except OSError as error: 
		print(error) 
		print("% s can not be created" % OutDir)

		
# Create new output folder for outputs when finding mid-Z
MidOutDir=OutDir+"/Z-mid"
if (exists(MidOutDir)==0):
	try:
		os.mkdir(MidOutDir)
		print("% s created successfully" % OutDir) 
	except OSError as error: 
		print(error) 
		print("% s can not be created" % OutDir)


if (HomeDir.endswith("/")==0):
	HomeDir=HomeDir+"/"
IJ.run("Close All")

#CC=ColorScale(cS,cE,colLen)
#colors=DrawScale(CC,600,100)
#out = os.path.join(OutDir,"ColorScale.tif")
#print("Saving to: "+out)
#IJ.saveAs(colors, "Tiff", out)
#IJ.saveAs("Results", os.path.join(OutDir,"Color Palette.csv"))
#IJ.selectWindow("Color Palette.csv");
#IJ.run("Close" )


listFiles = os.listdir(HomeDir)
# ct=0;
rt_main = ResultsTable()
# for (p=0; p < listPoints.length; p++){
#	title=listPoints[p];
for oib in listFiles:
	match1 = re.search("\.oib",oib)
	match2 = re.search("stack",oib)
	
	# if match1:
	# 	print("match1 passed!")
	# if match2:
	# 	print("match2 passed!")
	
	if match1 and match2 :
		print("Processing: " + HomeDir+oib)
		basename=re.sub(".oib","",oib)
		OutLineDir=OutDir+"/"+basename
		if (exists(OutLineDir)==0):
			try:
				os.mkdir(OutLineDir)
				print("% s created successfully" % OutDir) 
			except OSError as error: 
				print(error) 
				print("% s can not be created" % OutDir)
		
		IJ.run("Close All")
		IJ.run("Bio-Formats", "open="+HomeDir+oib+" color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		
		#sys.exit()		
		imp = IJ.getImage()
		#imp = IJ.getProcessor()
		# IJ.run(imp,"8-bit","")
		imp.show()
		#sys.exit()
		
		#Get Calibration
		cal = imp.getCalibration()
		print(type(cal))
		Xscale = cal.pixelWidth
		Yscale = cal.pixelHeight
		Unit = cal.unit
		
		#Set colors
		IJ.setForegroundColor(255, 255, 255)
		IJ.setBackgroundColor(0, 0, 0)
		
		#Record scales
		rt_main.incrementCounter()
		rt_main.addValue("Title", oib)
		rt_main.addValue("X Scale", Xscale)
		rt_main.addValue("Y Scale", Yscale)
		rt_main.addValue("Units", Unit)
		rt_main.show("Main Results")
		
		#get dimensions
		nS=imp.getNSlices()
		# print("number of slices:"+str(nS))
		rt_main.addValue("NSlices", nS)
		#Check if there are imcomplete scanning images using the TD channel
		ht=imp.height
		wd=imp.width
		# print(type(int(ht*0.9)))
		for i in range(nS):
			imp.setPosition(3,i+1,imp.getFrame())
			# imp.setSlice(i+1)
			imp.setRoi(int(wd*0.9),int(ht*0.9),int(wd*0.1),int(ht*0.1))
			proi = IJ.getImage().getRoi()
			stats = proi.getStatistics()
			if (stats.mean == 0):
				print("Mean of roi:"+str(stats.mean))
				# IJ.run(td_stack, "Draw", "slice")
				IJ.run(imp, "Delete Slice", "delete=slice")
		
		#Split channels
		realTtl = imp.getTitle()
		IJ.run(imp, "Split Channels", "")

		DAPI = "C1-"+realTtl
		EdU = "C2-"+realTtl
		TD ="C3-"+realTtl
		# print("title is "+title)
		
		#Get image from all channels and reverse the stacks
		IJ.selectWindow(DAPI)
		dapi_stack=IJ.getImage()
		dapi_stack.copyScale(imp)
		IJ.run(dapi_stack, "Blue", "")
		IJ.run(dapi_stack,"Reverse","")
		IJ.run(dapi_stack, "Divide...", "value=10 stack")
		
		IJ.selectWindow(EdU)
		edu_stack=IJ.getImage()
		IJ.run(edu_stack, "Red", "")
		IJ.run(edu_stack,"Reverse","")
		IJ.run(edu_stack, "Subtract...", "value=600 stack")
		IJ.run(edu_stack, "Multiply...", "value=1.8 stack")
		
		IJ.selectWindow(TD)
		td_stack=IJ.getImage()
		IJ.run(td_stack, "Grays", "")
		IJ.run(td_stack,"Reverse","")
		
				
		#Merge dapi and EdU channels
		dapi_avg= ZProjector.run(dapi_stack,"av")
		dapi_avg.copyScale(dapi_stack)
		dapi_avg.setTitle("dapi_sd")
		#IJ.run(dapi_avg,"8-bit","")
		IJ.run(dapi_avg,"Blue","")
		# dapi_avg.show()
		
		edu_max= ZProjector.run(edu_stack,"max")
		edu_max.copyScale(edu_stack)
		edu_max.setTitle("edu_max")
		#IJ.run(edu_max,"8-bit","")
		IJ.run(edu_max,"Red","")

		merged=RGBStackMerge.mergeChannels([edu_max,dapi_avg],False)
		IJ.run(merged, "RGB Color", "")
		merged.copyScale(imp)
		merged.show()
		# sys.exit()
		
		#Manual correction of image orientation
		degree = getAngle(merged,1)+90
		
		IJ.run(merged, "Rotate... ", "angle="+str(degree)+" grid=1 interpolation=Bilinear fill enlarge")
		IJ.run(edu_max, "Rotate... ", "angle="+str(degree)+" grid=1 interpolation=Bilinear fill enlarge")
		
		IJ.selectWindow(dapi_stack.getTitle())
		IJ.run(dapi_stack, "Rotate... ", "angle="+str(degree)+" grid=1 interpolation=Bilinear fill enlarge stack")
		
		IJ.selectWindow(edu_stack.getTitle())
		IJ.run(edu_stack, "Rotate... ", "angle="+str(degree)+" grid=1 interpolation=Bilinear fill enlarge stack")
		
		IJ.selectWindow(td_stack.getTitle())
		IJ.run(td_stack, "Rotate... ", "angle="+str(degree)+" grid=1 interpolation=Bilinear fill enlarge stack")
		
		# td_stack.show()
		# td_avg1.close()
		# td_avg.close()
		
		#Save files after rotation
		out = os.path.join(OutLineDir,basename+"-DAPI+EdU_XY.tif")
		print("Saving to: "+out)
		IJ.saveAs(merged, "Tiff", out)
		# merged.close()
		
		#out = os.path.join(OutLineDir,basename+"-maxEdU.tif")
		#print("Saving to: "+out)
		#IJ.run(edu_max,"Grays","")
		#IJ.saveAs(edu_max, "Tiff", out)
		

		td_min= ZProjector.run(td_stack,"min")
		td_min.copyScale(td_stack)
		td_min.setTitle("td_min")
		IJ.run(td_min,"8-bit","")
		out = os.path.join(OutLineDir,basename+"-TD-min.tif")
		print("Saving to: "+out)
		IJ.saveAs(td_min, "Tiff", out)
		td_min.close()
		td_stack.changes = False
		td_stack.close()
		
		#Stardist analysis using edu_max
		#Enlarge for better fitting
		edu_dup=edu_max.duplicate()
		edu_dup=edu_dup.resize(2*edu_max.width, 2*edu_max.height, "bilinear")
		edu_dup.setTitle("edu_dup")
		# IJ.run(edu_dup, "Subtract...", "value=25") #Remove background
		edu_dup.show()
		#Run stardist
		IJ.run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'edu_dup', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'95.0', 'percentileTop':'100.0', 'probThresh':'0.05', 'nmsThresh':'0.5', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]")
		
		#mark stardist outputs on merged
		IJ.selectWindow(merged.getTitle())
		IJ.run(merged, "Invert", "")
		IJ.setForegroundColor(255, 0, 0) #Select red as foreground color
		IJ.run("Line Width...", "line=1") #Line width

		rm=RoiManager.getRoiManager()
		nR=rm.getCount()
		for i in range(nR):
			roi=rm.getRoi(i)
			cen=roi.getContourCentroid()
			merged.getProcessor().drawDot(int(cen[0]/2),int(cen[1]/2))
		rm.close()
		#Save the results just by stardist
		out = os.path.join(StOutDir,basename+"-stardist.tif")
		print("Saving to: "+out)
		IJ.saveAs(merged, "Tiff", out)
		
		#Manual correction
		IJ.setForegroundColor(255, 255, 255)
		IJ.setBackgroundColor(255, 0, 0)
		IJ.run("Line Width...", "line=4")
		
		msgN="Correct stardist results using white!"
		IJ.setTool("polygon") 
		msg=WaitForUserDialog(msgN)
		msg.show()
		
		out = os.path.join(StOutDir,basename+"-manual.tif")
		print("Saving to: "+out)
		IJ.saveAs(merged, "Tiff", out)
		
		#Extract coordinates
		#Filter out red dots
		mask=ColorFilter(merged,[255,0,0])
		mask.show()
		IJ.run("Set Measurements...", "center bounding display redirect=None decimal=3");
		IJ.setRawThreshold(mask, 255, 255, "")
		IJ.run(mask, "Analyze Particles...", "pixel show=Nothing display clear")
		
		out = os.path.join(CoOutDir,basename+".csv")
		print("Saving to: "+out)
		IJ.saveAs("Results", out)


		# sys.exit()
				
		#Color-code the EdU stack
		edu_cc=ColorCode(imp=edu_stack,start=1,end=edu_stack.getNSlices(),cS=blue,cE=orange,limit=1000)
		IJ.run(edu_stack, "Red", "") #Recover the red LUT
		edu_cc.setTitle("edu_cc")
		edu_cc.copyScale(edu_stack)
		edu_cc.show()

		#sys.exit()
		# continue
		#Make maximum projections at XY direction
		edu_cc_proj = ZProjector.run(edu_cc,"max")
		edu_cc_proj.copyScale(edu_stack)
		edu_cc_proj.setTitle("edu_cc_proj")
		#Save and close
		out = os.path.join(OutLineDir,basename+"-EdU-XY.tif")
		edu_cc_proj.show()
		print("Saving EdU CC XY proj to: "+out)
		IJ.saveAs(edu_cc_proj, "Tiff", out)
		edu_cc_proj.close()
		
		#Create XY mask from dapi channel
		dapi_proj_XY=ZProjector.run(dapi_stack,"avg")
		dapi_mask_XY=GetMask(dapi_proj_XY)
		dapi_mask_XY.setTitle("dapi_mask_XY")
		dapi_mask_XY.show()
		#Save and close
		out = os.path.join(OutLineDir,basename+"-XY_mask.tif")
		print("Saving XY mask to: "+out)
		IJ.saveAs(dapi_mask_XY, "Tiff", out)
		#sys.exit()
		dapi_mask_XY.close()
		
		#Reslice the image to get XZ projections
		IJ.run(dapi_stack, "Reslice [/]...", "output=5.000 start=Top")
		dapi_stack_XZ=IJ.getImage()
		dapi_stack_XZ.setTitle("dapi_stack_XZ")
		#Maximum projection
		dapi_proj_XZ=ZProjector.run(dapi_stack_XZ,"avg")
		dapi_proj_XZ.setTitle("dapi_proj_XZ")
		# Add some distance to edges, in order to show the boundary in z direction
		IJ.run(dapi_proj_XZ, "Canvas Size...", "width="+str(dapi_proj_XZ.width)+" height="+str(dapi_proj_XZ.height+30)+" position=Center zero")
		dapi_proj_XZ.show()
		
		#Binary XZ mask from DAPI channel
		dapi_mask_XZ=GetMask(dapi_proj_XZ)
		dapi_mask_XZ.setTitle("dapi_mask_XZ")
		dapi_mask_XZ.show()
		#Save and close
		out = os.path.join(OutLineDir,basename+"-XZ_mask.tif")
		print("Saving XZ mask to: "+out)
		IJ.saveAs(dapi_mask_XZ, "Tiff", out)
		dapi_mask_XZ.close()
		
		#Repeat reslicing for original EdU
		IJ.run(edu_stack, "Reslice [/]...", "output=5.000 start=Top")
		edu_stack_XZ=IJ.getImage()
		edu_stack_XZ.setTitle("edu_stack_XZ")
		#Maximum projection
		edu_proj_XZ=ZProjector.run(edu_stack_XZ,"max")
		edu_proj_XZ.setTitle("edu_proj_XZ")
		# Add some distance to edges, in order to show the boundary in z direction
		IJ.run(edu_proj_XZ, "Canvas Size...", "width="+str(edu_proj_XZ.width)+" height="+str(edu_proj_XZ.height+30)+" position=Center zero")
		
		#DAPI+EdU in XZ direction
		merged_XZ=RGBStackMerge.mergeChannels([edu_proj_XZ,dapi_proj_XZ],True)
		IJ.run(merged_XZ, "RGB Color", "")
		merged_XZ.copyScale(dapi_stack_XZ)
		#Save and close
		out = os.path.join(OutLineDir,basename+"-DAPI+EdU_XZ.tif")
		print("Saving merged XZ to: "+out)
		IJ.saveAs(merged_XZ, "Tiff", out)
		merged_XZ.close()
		dapi_proj_XZ.close()
		edu_proj_XZ.close()
		#sys.exit()
		
		#Repeat reslicing for color-coded EdU
		IJ.run(edu_cc, "Reslice [/]...", "output=5.000 start=Top")
		edu_cc_XZ=IJ.getImage()
		edu_cc_XZ.setTitle("edu_cc_XZ")
		#Maximum projection
		edu_cc_proj_XZ=ZProjector.run(edu_cc_XZ,"max")
		edu_cc_proj_XZ.setTitle("edu_cc_proj_XZ")
		# Add some distance to edges, in order to show the boundary in z direction
		IJ.run(edu_cc_proj_XZ, "Canvas Size...", "width="+str(edu_cc_proj_XZ.width)+" height="+str(edu_cc_proj_XZ.height+30)+" position=Center zero")
		edu_cc_proj_XZ.show()

		#Saving and closing XZ directions
		out = os.path.join(OutLineDir,basename+"-EdU_XZ.tif")
		print("Saving to: "+out)
		IJ.saveAs(edu_cc_proj_XZ, "Tiff", out)
		edu_cc_proj_XZ.close()
	
		edu_cc_XZ.close()
		dapi_stack_XZ.close()
		edu_stack_XZ.close()
		

		#Reslice the image to get YZ projections
		IJ.run(dapi_stack, "Reslice [/]...", "output=5.000 start=Left rotate")
		dapi_stack_YZ=IJ.getImage()
		dapi_stack_YZ.setTitle("dapi_stack_YZ")
		#Maximum projection
		dapi_proj_YZ=ZProjector.run(dapi_stack_YZ,"avg")
		dapi_proj_YZ.setTitle("dapi_proj_YZ")
		# Add some distance to edges, in order to show the boundary in z direction
		IJ.run(dapi_proj_YZ, "Canvas Size...", "width="+str(dapi_proj_YZ.width+30)+" height="+str(dapi_proj_YZ.height)+" position=Center zero")
		
		#Binary mask from DAPI channel
		dapi_mask_YZ=GetMask(dapi_proj_YZ)
		dapi_mask_YZ.setTitle("dapi_mask_YZ")
		#Save and close
		out = os.path.join(OutLineDir,basename+"-YZ_mask.tif")
		print("Saving YZ mask to: "+out)
		IJ.saveAs(dapi_mask_YZ, "Tiff", out)
		dapi_mask_XZ.close()
		
		
		#Repeat reslicing for original EdU
		IJ.run(edu_stack, "Reslice [/]...", "output=5.000 start=Left rotate")
		edu_stack_YZ=IJ.getImage()
		edu_stack_YZ.setTitle("edu_stack_YZ")
		#Maximum projection
		edu_proj_YZ=ZProjector.run(edu_stack_YZ,"max")
		edu_proj_YZ.setTitle("edu_proj_YZ")
		# Add some distance to edges, in order to show the boundary in z direction
		IJ.run(edu_proj_YZ, "Canvas Size...", "width="+str(edu_proj_YZ.width+30)+" height="+str(edu_proj_YZ.height)+" position=Center zero")
		
		#DAPI+EdU in YZ direction
		merged_YZ=RGBStackMerge.mergeChannels([edu_proj_YZ,dapi_proj_YZ],True)
		IJ.run(merged_YZ, "RGB Color", "")
		merged_YZ.copyScale(dapi_stack_YZ)
		merged_YZ.show()
		#Save and close
		out = os.path.join(OutLineDir,basename+"-DAPI+EdU_YZ.tif")
		print("Saving to: "+out)
		IJ.saveAs(merged_YZ, "Tiff", out)
		merged_YZ.close()
		edu_proj_YZ.close()
		dapi_proj_YZ.close()

		#Repeat reslicing for color-coded EdU
		IJ.run(edu_cc, "Reslice [/]...", "output=5.000 start=Left rotate")
		edu_cc_YZ=IJ.getImage()
		edu_cc_YZ.setTitle("edu_cc_YZ")
		#Maximum projection
		edu_cc_proj_YZ=ZProjector.run(edu_cc_YZ,"max")
		edu_cc_proj_YZ.setTitle("edu_cc_proj_YZ")
		# Add some distance to edges, in order to show the boundary in z direction
		IJ.run(edu_cc_proj_YZ, "Canvas Size...", "width="+str(edu_cc_proj_YZ.width+30)+" height="+str(edu_cc_proj_YZ.height)+" position=Center zero")
		edu_cc_proj_YZ.show()

		#Saving and closing YZ directions
		out = os.path.join(OutLineDir,basename+"-EdU_YZ.tif")
		print("Saving to: "+out)
		IJ.saveAs(edu_cc_proj_YZ, "Tiff", out)
		edu_cc_proj_YZ.close()
		edu_cc_YZ.close()
		dapi_stack_YZ.close()
		edu_stack_YZ.close()

		print("Finished for the current image stack.")
		
		IJ.run("Close All")	
		# pass
	# matches("test.oib1", ".*oib.*");
out = os.path.join(OutDir,"Statistics.csv")
print("Saving to: "+out)
rt_main.save(out)
print("All Finished!")
