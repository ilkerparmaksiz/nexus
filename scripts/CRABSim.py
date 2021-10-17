import multiprocessing as mtp
import numpy as np
import matplotlib.pyplot as plt
import h5py as h
import time
import pickle
import os
import random
#Create a custom 2d Histogram
saveFile="/media/ilker/writable/DATA/pic/"
def Hist2d(title,Energys,Times,binss,xlimit,ylimit):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.hist2d(Energys,Times,bins=binss,cmap=plt.cm.BuGn_r)
    plt.title(title)

    plt.ylabel("Energy keV",fontsize=26)
    plt.xlabel("Time us",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlim(xlimit)
    plt.ylim(ylimit)
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(saveFile+title + '.png')
    plt.show()

#Create a custom 1d Histogram
def Hist1d(title,Energys,Bins,xlimit,ylimit,FileSave,limits=False):
    fig, ax = plt.subplots(figsize=(8,8))
    plt.hist(Energys,bins=Bins,alpha=0.7,label=title)
    plt.ylabel("Counts",fontsize=26)
    plt.xlabel("Energy keV",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='upper right',fontsize=24,shadow=True, fancybox=True)
    #plt.semilogy()
    if(limits):
        plt.xlim(xlimit)
        plt.ylim(ylimit)
    plt.grid(True)
    plt.tight_layout()
    if(FileSave):
        savetit=saveFile+title + '.png'
        print(f"saving the picture at {savetit}")
        plt.savefig(savetit)
    else:
        plt.show()

#Read the HDF5 file and return the dataset
def ReadFile(file):
    data = h.File(file,'r')
    return data

#Create a custom 2d Plot
def Plot2d(title,x,y):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.plot(x,y)
    plt.title(title)
    plt.ylabel("Z mm",fontsize=26)
    plt.xlabel("Y mm",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

#This draws two circles on the top of an other
def Circle(pltX,pltY,title,FileSave,LarCirR=75,SmallCirR=35,xlimit=100,ylimit=100):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.title(title)
    plt.plot(pltX,pltY)
    plt.scatter(pltX,pltY)
    circle1 = plt.Circle((0, 0), SmallCirR, color='r',fill=False)
    circle2 = plt.Circle((0, 0), LarCirR, color='k',fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    plt.xlim(-xlimit,xlimit)
    plt.ylim(-ylimit,ylimit)

    plt.ylabel("Z (mm)",fontsize=26)
    plt.xlabel("Y (mm)",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.grid(True)
    plt.tight_layout()
    if(FileSave):
        savetit=saveFile+title + '.png'
        print(f"saving the picture at {savetit}")
        plt.savefig(savetit)
    else:
        plt.show()

#This function runs and collects data if processes > 1 it will do multiprocessing
def RunEvents(file,Processes=4):
    data=ReadFile(file)
    Energys=[]
    TotalEvents=int(data['MC']['configuration'][2][1])
    print(f"Currently there are {TotalEvents}")
    print(f"Processing Your Events at {file}")
    if(Processes>1):
        XEvent=int(TotalEvents/Processes)
        p=[]
        E=mtp.Array('d',range(TotalEvents))
        for i in range(0,Processes):
            p.append(mtp.Process(target=getEnergys,args=(data,i*XEvent,(i+1)*XEvent,E)))
        for x in p:
            x.start()
        for k in p:
            k.join()

        Energys=np.array(E)

        print(f"There are {len(Energys)} proccessed")
        return Energys
    else:
        Energys=np.array(getEnergysSingle(data,0,TotalEvents))
    return Energys

#Collects only Energy Spectrum of the specified events
def getEnergys(data,xRLower,xRHigher,E):
    for x in range(xRLower,xRHigher):
        Current_Event = x
        Current_Hit_Mask = data["MC"]['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event

        Electron_Mask = data["MC"]['particles'][Current_Particle_Mask]['particle_name'] == b'e-'

        Electron_PIDS = data["MC"]['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']

        A = data["MC"]['hits'][Current_Hit_Mask]['particle_id']
        Hit_Electron_Maks = np.in1d(A, Electron_PIDS)

        E[x]=(data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]['energy'].sum()*1e3)
#Collects Energy for non multiprocessing
def getEnergysSingle(data,xRLower,xRHigher):
    Energys=[]
    for x in range(xRLower,xRHigher):
        Current_Event = x
        Current_Hit_Mask = data["MC"]['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event

        Electron_Mask = data["MC"]['particles'][Current_Particle_Mask]['particle_name'] == b'e-'

        Electron_PIDS = data["MC"]['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']

        A = data["MC"]['hits'][Current_Hit_Mask]['particle_id']
        Hit_Electron_Maks = np.in1d(A, Electron_PIDS)

        Energys.append(data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]['energy'].sum()*1e3)
    return Energys

def getTracks(file,Current_Event=0):
    data = h.File(file,'r')

    Current_Hit_Mask = data["MC"]['hits']['event_id'] == Current_Event

    Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event

    Electron_Mask = data["MC"]['particles'][Current_Particle_Mask]['particle_name'] == b'e-'

    Electron_PIDS = data["MC"]['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']

    A = data["MC"]['hits'][Current_Hit_Mask]['particle_id']
    Hit_Electron_Maks = np.in1d(A, Electron_PIDS)
    ElectronData=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]
    PIDWithCounts={}

    for i in ElectronData["particle_id"]:
        PID_mask=ElectronData["particle_id"]==i
        count = np.count_nonzero(PID_mask)
        PIDWithCounts[count]=PID_mask

    max_key=max(list(PIDWithCounts.keys()))
    ElectronsWithMaxTrack=ElectronData[PIDWithCounts[max_key]]
    Tracks=np.array([ElectronsWithMaxTrack["x"],ElectronsWithMaxTrack["y"],ElectronsWithMaxTrack["z"]])

    return Tracks


#Gets the Tracks and as well as energies
def RunEventsDic(file,npFileName,Processes=4,TotalLimit=0):
    data=ReadFile(file)
    theEvents={}
    SimEvents=int(data['MC']['configuration'][1][1])
    if TotalLimit==0:
        TotalEvents=int(data['MC']['configuration'][2][1])
    else:
        TotalEvents=TotalLimit
    print(f"Processing Your Events at {file}")
    print(f"Number of SimEvents is {SimEvents}")
    print(f"Number Interacted Events is {TotalEvents} events")

    if(Processes>1):
        manager=mtp.Manager()
        theEvents=manager.dict()

        XEvent=int(TotalEvents/Processes)
        p=[]
        for i in range(0,Processes):
            p.append(mtp.Process(target=EnergyAndTracks,args=(data,i*XEvent,(i+1)*XEvent,theEvents)))
        for x in p:
            x.start()
        for k in p:
            k.join()


        print(f"There are {len(theEvents)} proccessed")

    else:
        EnergyAndTracks(data,0,TotalEvents,theEvents)

    npFileName=npFileName+".data"
    #SaveFile(npFileName,theEvents)

    with open(npFileName, 'wb') as outfile:
        pickle.dump(theEvents, outfile, protocol=pickle.HIGHEST_PROTOCOL)
    return theEvents



def EnergyAndTracks(data,xRLower,xRHigher,theEvents,EventLimit=0,fudicalR=35,fdcount=True,AllEvents=False):
    TotalEvents=xRHigher
    if(EventLimit>0):
        TotalEvents=EventLimit


    FudicalCount=0
    for Current_Event in range(xRLower,TotalEvents):
        Current_Hit_Mask = data['MC']['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event
        Hits_PIDs=data['MC']['hits'][Current_Hit_Mask]['particle_id']
        Current_Particles=data['MC']['particles'][Current_Particle_Mask]

        Electron_Mask = Current_Particles['particle_name'] == b'e-'
        Gamma_Mask = Current_Particles['particle_name'] == b'gamma'


        Electron_MIDS = Current_Particles[Electron_Mask]['mother_id']
        Gamma_PIDS = Current_Particles[Gamma_Mask]['particle_id']


        Mothers=np.in1d(Gamma_PIDS,Electron_MIDS)

        MotherEnergy=Current_Particles[Gamma_Mask][Mothers]['kin_energy']*1e3


        Electron_PIDS = data['MC']['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']


        Hit_Electron_Maks = np.in1d(Hits_PIDs, Electron_PIDS)
        ElectronData=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]

        TotalEventEnergy=ElectronData["energy"].sum()*1e3

        Xhits=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]["x"]
        Yhits=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]["y"]
        Zhits=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]["z"]

        #Checking the Radius if in the boundaries
        fRadius=np.sqrt(Yhits*Yhits+Zhits*Zhits)
        fRadius_Mask=fRadius<=fudicalR
        tracklength=len(Yhits)
        TrksCntIntheFudi=np.count_nonzero(fRadius_Mask)
        RatioTrackLength=(TrksCntIntheFudi/tracklength)*100
        if(not AllEvents):
            if(TotalEventEnergy<=300):
                if(RatioTrackLength<=95) :
                    continue
        FudicalCount+=1

        Tracks=np.array([Xhits,Yhits,Zhits])
        theEvents[Current_Event]=[Tracks,TotalEventEnergy,TrksCntIntheFudi,MotherEnergy]


    if(fdcount):
        print(f" There are only {FudicalCount} out of {TotalEvents} in the fudical volume")


def PlotRandomEvents(theEvents,NPlots,LowE,HighE,pretitle,FileSave=True):

    QualfEnergys=[]
    for Event in theEvents:
        if(theEvents[Event][1]>=LowE and theEvents[Event][1]<=HighE):
            QualfEnergys.append(Event)
    print(f"Events at ({LowE} - {HighE}) -> {len(QualfEnergys)}")
    for plot in range(0,NPlots):
        Event=random.choice(QualfEnergys)
        title=str(pretitle) + "_Tracks_" + str(round(theEvents[Event][1],2)) + "keV"
        Circle(theEvents[Event][0][1],theEvents[Event][0][2],title,FileSave)


def PlotEnergySpec(theEvents,title,binss=np.arange(1,1000,10),FileSave=True):
    Energys=[]
    for Event in theEvents:
        Energys.append(theEvents[Event][1])
    Hist1d(title,Energys,binss,0,0,FileSave)

def PlotMotherEnergySpec(theEvents,title,binss=np.arange(1,1000,10),FileSave=True):
    Energys2=[]
    for Event in theEvents:
        MotherE=theEvents[Event][3]
        Energys2.extend(MotherE)
    Hist1d(title,Energys2,binss,0,0,FileSave)

def LoadPickleFile(file):
    with open(file, 'rb') as f:
        theEvents=pickle.load(f)
    return theEvents

def Hist2dv2(title,Energys,Times,binss,xlimit,ylimit):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.hist2d(Energys,Times,bins=binss,cmap=plt.cm.BuGn_r)
    plt.title(title)

    plt.ylabel("Energy keV",fontsize=26)
    plt.xlabel("Time us",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.xlim(xlimit)
    plt.ylim(ylimit)
    plt.grid(True)
    plt.tight_layout()

    plt.savefig('/home/ilker/Dropbox/nexus/build/source/'+title + '.png')
    plt.show()

def Hist1dv2(title,Energys,Bins,xlabel,ylabel):
    fig, ax = plt.subplots(figsize=(8,8))
    plt.hist(Energys,bins=Bins,alpha=0.7,label=title)
    plt.ylabel(ylabel,fontsize=26)
    plt.xlabel(xlabel,fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='upper right',fontsize=24,shadow=True, fancybox=True)
    #plt.semilogy()

    plt.grid(True)
    plt.tight_layout()

    #plt.savefig('/Users/austinmcdonald/Desktop/'+title + '.png')
    plt.show()

def Plot2dv2(title,x,y):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.plot(x,y)
    plt.scatter(x,y)
    plt.title(title)
    plt.ylabel("Z mm",fontsize=26)
    plt.xlabel("Y mm",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def Plot3dv2(title,x,y,z):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, 'gray')
    ax.scatter3D(x, y, z, c=z, cmap='Greens');
    plt.title(title)
    plt.ylabel("Y mm",fontsize=26)
    plt.xlabel("X mm",fontsize=26)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def Circlev2(pltX,pltY,title,LarCirR=75,SmallCirR=35,xlimit=100,ylimit=100,FileSave=False):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.title(title)
    plt.plot(pltX,pltY)
    plt.scatter(pltX,pltY)
    circle1 = plt.Circle((0, 0), SmallCirR, color='r',fill=False)
    circle2 = plt.Circle((0, 0), LarCirR, color='k',fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    plt.xlim(-xlimit,xlimit)
    plt.ylim(-ylimit,ylimit)

    plt.ylabel("Z (mm)",fontsize=26)
    plt.xlabel("Y (mm)",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.grid(True)
    plt.tight_layout()
    if(FileSave):
        plt.savefig('/home/ilker/Pictures/Sim/'+title + '.png')
    plt.show()


def EnergyAndTracksorig(data,EventLimit=0,fudicalR=35,fdcount=True,AllEvents=False):
    #data = h.File(file,'r')
    TotalEvents=int(data['MC']['configuration'][2][1])
    if(EventLimit>0):
        TotalEvents=EventLimit

    theEvents={}

    FudicalCount=0
    ElectronHits={}


    for Current_Event in range(0,TotalEvents):

        Current_Hit_Mask = data['MC']['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event
        Hits_PIDs=data['MC']['hits'][Current_Hit_Mask]['particle_id']
        Current_Particles=data['MC']['particles'][Current_Particle_Mask]

        Electron_Mask = Current_Particles['particle_name'] == b'e-'
        #Gamma_Mask = Current_Particles['particle_name'] == b'gamma'


        Electron_MIDS = Current_Particles[Electron_Mask]['mother_id']
        #Gamma_PIDS = Current_Particles[Gamma_Mask]['particle_id']


        #Mothers=np.in1d(Gamma_PIDS,Electron_MIDS)

        #MotherEnergy=Current_Particles[Gamma_Mask][Mothers]['kin_energy']
        '''print (f"Information about GammaRays for event {Current_Event}")
        
        print ("All the Gammas in the particles")
        print(data['MC']['particles'][Current_Particle_Mask][Gamma_Mask]['kin_energy']*1e3)
        
        print("Gammas Mother to Electrons in the event")
        print(MotherEnergy*1000)
        '''

        Electron_PIDS = data['MC']['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']


        Hit_Electron_Maks = np.in1d(Hits_PIDs, Electron_PIDS)
        ElectronData=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]

        TotalEventEnergy=ElectronData["energy"].sum()*1e3
        #print (f"Total Energy for event {Current_Event} is {TotalEventEnergy} ")

        Xhits=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]["x"]
        Yhits=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]["y"]
        Zhits=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]["z"]

        #Checking the Radius if in the boundaries
        fRadius=np.sqrt(Yhits*Yhits+Zhits*Zhits)
        fRadius_Mask=fRadius<=fudicalR
        tracklength=len(Yhits)
        TrksCntIntheFudi=np.count_nonzero(fRadius_Mask)
        RatioTrackLength=(TrksCntIntheFudi/tracklength)*100
        if(not AllEvents):
            if(RatioTrackLength<=95):
                continue
        FudicalCount=FudicalCount+1

        Tracks=np.array([Xhits,Yhits,Zhits])
        theEvents[Current_Event]=[Tracks,TotalEventEnergy,TrksCntIntheFudi]

    if(fdcount):
        print(f" There are only {FudicalCount} out of {TotalEvents} in the fudical volume")

    return theEvents

def PlotRandomEvents(theEvents,NPlots,LowE,HighE,EventLimit=0,FileSave=False):

    QualfEnergys=[]
    for Event in theEvents:
        if(theEvents[Event][1]>=LowE and theEvents[Event][1]<=HighE):
            QualfEnergys.append(Event)
    print(len(QualfEnergys))
    for plot in range(0,NPlots):
        Event=random.choice(QualfEnergys)
        title="Track with " + str(theEvents[Event][1]) + " keV"
        Circle(theEvents[Event][0][1],theEvents[Event][0][2],title,FileSave)



def PlotEnergySpec(theEvents,title,binss=np.arange(1,500,10)):
    Energys=[]
    for Event in theEvents:
        Energys.append(theEvents[Event][1])

    print (len(Energys))
    Hist1d(title,Energys,binss,"Counts","Energy (keV)",True)

def EnergyAndTracksv2(file,particle=b'e-',EventLimit=0):
    data = h.File(file,'r')
    TotalEvents=int(data['MC']['configuration'][1][1])
    IntEvnts=int(data['MC']['configuration'][2][1])
    RatioInterTotal=IntEvnts/TotalEvents

    if(EventLimit>0):
        TotalEvents=EventLimit

    theEvents={}
    print(f"Total Events for {particle} is {TotalEvents}")
    print (f"Interacted Events for {particle} is {IntEvnts}")
    print (f"Ratio Interc/Total is {RatioInterTotal}")

    FudicalCount=0

    for Current_Event in range(0,IntEvnts):

        Current_Hit_Mask = data['MC']['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event
        Hits_PIDs=data['MC']['hits'][Current_Hit_Mask]['particle_id']
        Current_Particles=data['MC']['particles'][Current_Particle_Mask]

        Electron_Mask = Current_Particles['particle_name'] == particle


        Electron_MIDS = Current_Particles[Electron_Mask]['mother_id']


        Electron_PIDS = data['MC']['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']


        Hit_Electron_Maks = np.in1d(Hits_PIDs, Electron_PIDS)
        ElectronData=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]

        TotalEventEnergy=ElectronData["energy"].sum()*1e3
        PartID=np.unique(ElectronData["particle_id"])
        ElectronHits={}
        if (len(PartID)==0):
            continue
        for ID in PartID:
            ElectronTrackMask=ElectronData["particle_id"]==ID
            Xhits=ElectronData[ElectronTrackMask]["x"]
            Yhits=ElectronData[ElectronTrackMask]["y"]
            Zhits=ElectronData[ElectronTrackMask]["z"]
            T_Xhit=np.diff(Xhits).sum()
            T_Yhit=np.diff(Yhits).sum()
            T_Zhit=np.diff(Zhits).sum()
            length=math.sqrt(T_Xhit*T_Xhit+T_Yhit*T_Yhit+T_Zhit*T_Zhit)
            if(length==0):
                continue
            Tracks=np.array([Xhits,Yhits,Zhits])
            ElectronHits[ID]=[Tracks,length]

        theEvents[Current_Event]=[ElectronHits,TotalEventEnergy]

    return theEvents

def EnergyAndTrackscf137(file,mother=b'gamma',EventLimit=0):
    data = h.File(file,'r')

    TotalEvents=int(data['MC']['configuration'][1][1])
    IntEvnts=int(data['MC']['configuration'][2][1])
    RatioInterTotal=IntEvnts/TotalEvents

    if(EventLimit>0):
        TotalEvents=EventLimit

    print (f"Total Events for Cs137 is {TotalEvents}")
    print (f"Interacted Events for Cs137 is {IntEvnts}")
    print (f"Ratio Interc/Total is {RatioInterTotal}")

    theEvents={}

    FudicalCount=0

    for Current_Event in range(0,IntEvnts):

        Current_Hit_Mask = data['MC']['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event
        Hits_PIDs=data['MC']['hits'][Current_Hit_Mask]['particle_id']
        Current_Particles=data['MC']['particles'][Current_Particle_Mask]

        Electron_Mask = Current_Particles['particle_name'] == b'e-'
        Electron_MIDS = Current_Particles[Electron_Mask]['mother_id']

        Gamma_Mask = Current_Particles['particle_name'] == mother
        Gamma_PrimMother_Mask = Current_Particles[Gamma_Mask]['mother_id']==2
        Gamma_PIDS = Current_Particles[Gamma_Mask][Gamma_PrimMother_Mask]['particle_id']


        MothersMask=np.in1d(Electron_MIDS,Gamma_PIDS)
        Electron_PIDS = data['MC']['particles'][Current_Particle_Mask][Electron_Mask][MothersMask]['particle_id']

        Hit_Electron_Maks = np.in1d(Hits_PIDs, Electron_PIDS)
        ElectronData=data['MC']['hits'][Current_Hit_Mask][Hit_Electron_Maks]

        TotalEventEnergy=ElectronData["energy"].sum()*1e3
        PartID=np.unique(ElectronData["particle_id"])
        ElectronHits={}

        if (len(PartID)==0):
            continue

        for ID in PartID:
            ElectronTrackMask=ElectronData["particle_id"]==ID

            Xhits=ElectronData[ElectronTrackMask]["x"]
            Yhits=ElectronData[ElectronTrackMask]["y"]
            Zhits=ElectronData[ElectronTrackMask]["z"]
            T_Xhit=np.diff(Xhits).sum()
            T_Yhit=np.diff(Yhits).sum()
            T_Zhit=np.diff(Zhits).sum()
            length=math.sqrt(T_Xhit*T_Xhit+T_Yhit*T_Yhit+T_Zhit*T_Zhit)
            DistanceToFirstPoint=math.sqrt(Xhits[0]*Xhits[0]+Yhits[0]*Yhits[0]+Zhits[0]*Zhits[0])

            if(length==0):
                continue
            Tracks=np.array([Xhits,Yhits,Zhits])
            ElectronHits[ID]=[Tracks,length,DistanceToFirstPoint]

        theEvents[Current_Event]=[ElectronHits,TotalEventEnergy,DistanceToFirstPoint]

    return theEvents

def Sr90EandTrk(file,EventLimit=0):
    data = h.File(file,'r')
    TotalEvents=int(data['MC']['configuration'][1][1])
    IntEvnts=int(data['MC']['configuration'][2][1])
    RatioInterTotal=IntEvnts/TotalEvents
    print (f"Total Events for Sr90 is {TotalEvents}")
    print (f"Interacted Events for Sr90 is {IntEvnts}")
    print (f"Ratio Interc/Total is {RatioInterTotal}")
    if(EventLimit>0):
        TotalEvents=EventLimit

    theEvents={}

    FudicalCount=0

    for Current_Event in range(0,IntEvnts):
        Current_Hit_Mask = data['MC']['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event
        Hits_PIDs=data['MC']['hits'][Current_Hit_Mask]['particle_id']
        Current_Particles=data['MC']['particles'][Current_Particle_Mask]

        Electron_Mask = Current_Particles['particle_name'] == b'e-'
        Electron_MIDS = Current_Particles[Electron_Mask]['mother_id']
        Electron_PIDS = data['MC']['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']

        Sr90Beta_Mask = Electron_MIDS == 1
        Y90Beta_Mask = Electron_MIDS == 2

        Sr90BetaIDs=Current_Particles[Electron_Mask][Sr90Beta_Mask]['particle_id']
        Y90BetaIDs=Current_Particles[Electron_Mask][Y90Beta_Mask]['particle_id']

        Hit_Sr90Beta_Mask = np.in1d(Hits_PIDs, Sr90BetaIDs)
        Hit_Y90Beta_Mask = np.in1d(Hits_PIDs, Y90BetaIDs)
        Sr90Betas=data['MC']['hits'][Current_Hit_Mask][Hit_Sr90Beta_Mask]
        Y90Betas=data['MC']['hits'][Current_Hit_Mask][Hit_Y90Beta_Mask]


        Sr90E=Sr90Betas['energy'].sum()*1000
        Y90E=Y90Betas['energy'].sum()*1000



        Sr90BetaHits={}
        Y90BetaHits={}

        # BetaHits from Sr90
        Sr90Xhits=Sr90Betas["x"]
        Sr90Yhits=Sr90Betas["y"]
        Sr90Zhits=Sr90Betas["z"]
        Sr90T_Xhit=np.diff(Sr90Xhits).sum()
        Sr90T_Yhit=np.diff(Sr90Yhits).sum()
        Sr90T_Zhit=np.diff(Sr90Zhits).sum()
        Sr90DistancetoFirstPoint=math.sqrt(Sr90Xhits[0]*Sr90Xhits[0]+Sr90Yhits[0]*Sr90Yhits[0]+Sr90Zhits[0]*Sr90Zhits[0])
        Sr90length=math.sqrt(Sr90T_Xhit*Sr90T_Xhit+Sr90T_Yhit*Sr90T_Yhit+Sr90T_Zhit*Sr90T_Zhit)
        Sr90Tracks=np.array([Sr90Xhits,Sr90Yhits,Sr90Zhits])
        Sr90BetaHits[Current_Event]=[Sr90Tracks,Sr90length,Sr90E,Sr90DistancetoFirstPoint]


        #BetaHits from Y90
        Y90Xhits=Y90Betas["x"]
        Y90Yhits=Y90Betas["y"]
        Y90Zhits=Y90Betas["z"]
        Y90T_Xhit=np.diff(Y90Xhits).sum()
        Y90T_Yhit=np.diff(Y90Yhits).sum()
        Y90T_Zhit=np.diff(Y90Zhits).sum()
        Y90length=math.sqrt(Y90T_Xhit*Y90T_Xhit+Y90T_Yhit*Y90T_Yhit+Y90T_Zhit*Y90T_Zhit)
        Y90Tracks=np.array([Y90Xhits,Y90Yhits,Y90Zhits])
        Y90BetaHits[Current_Event]=[Y90Tracks,Y90length,Y90E]


        theEvents[Current_Event]=[Sr90BetaHits,Y90BetaHits]

    return theEvents


def PlotLongTracks(theEvents,FileSave=False):

    count=0
    for Event in theEvents:
        ElcTracks=theEvents[Event][0]
        title="Track with " + str(theEvents[Event][1]) + " keV"
        for ID in ElcTracks:
            if(ElcTracks[ID][1]>30 and count<5):
                #Plot3d(title,ElcTracks[ID][0][0],ElcTracks[ID][0][1],ElcTracks[ID][0][2])
                Circle(ElcTracks[ID][0][1],ElcTracks[ID][0][2],title,FileSave)
                count=count+1


def PlotLength(theEvents,title,binss=np.arange(1,500,10)):
    Distance=[]
    for Event in theEvents:
        ElcTracks=theEvents[Event][0]
        for ID in ElcTracks:
            Distance.append(ElcTracks[ID][1])
    Hist1d(title,Distance,binss,"Counts","Track Length (mm)")


def PlotsforSr90(theEvents,title,lenbins=np.arange(1,80,10),ebins=np.arange(1,3000,60),enrg=1000,trckthreshold=30,pltlimit=2,FileSave=False):
    Sr90Energys=[]
    Y90Energys=[]
    Sr90Distances=[]
    Sr90FirstDistances=[]
    Y90Distances=[]

    Sr90cnt=0
    Y90cnt=0
    print (f"Number of Interacted Events is {len(theEvents)}")
    for Event in theEvents:
        Sr90Beta=theEvents[Event][0]
        Y90Beta=theEvents[Event][1]
        for PID in Sr90Beta:
            Sr90Energys.append(Sr90Beta[PID][2])

            Sr90Distances.append(Sr90Beta[PID][1])
            Sr90FirstDistances.append(Sr90Beta[PID][3])
            title="Track with Sr90 " + str(Sr90Beta[PID][2]) + " keV"
            if(Sr90Beta[PID][1]>trckthreshold and Sr90cnt<=pltlimit):
                Circle(Sr90Beta[PID][0][1],Sr90Beta[PID][0][2],title,FileSave)
                Sr90cnt=Sr90cnt+1

        for YPID in Y90Beta:
            Y90Energys.append(Y90Beta[YPID][2])
            if(Y90Beta[YPID][2]>enrg):
                Y90Distances.append(Y90Beta[YPID][1])
                title="Track with Y90 " + str(Y90Beta[YPID][2]) + " keV"
                if(Y90Beta[YPID][1]>trckthreshold):
                    if(Y90cnt<=pltlimit):
                        Circle(Y90Beta[YPID][0][1],Y90Beta[YPID][0][2],title,FileSave)
                    Y90cnt=Y90cnt+1
    print(f"Y90 event with energy above {enrg} keV is {len(Y90Distances)} ")
    print(f"Number of Events with {trckthreshold} mm and higher tracks is {Y90cnt}")
    Hist1d("Sr90 Energy Spectrum",Sr90Energys,ebins,"Energy (keV)","Counts")
    Hist1d("Y90 Energy Spectrum",Y90Energys,ebins,"Energy (keV)","Counts")
    BEnergys=[]
    BEnergys.extend(Y90Energys)
    BEnergys.extend(Sr90Energys)
    Hist1d("S90 and Y90 Energy Spectrum Combined",BEnergys,ebins,"Energy (keV)","Counts")
    Hist1d("Sr90 First Hit Distance",Sr90FirstDistances,np.arange(1,20,2),"Length (mm)","Counts")
    Hist1d("Sr90Beta Track",Sr90Distances,lenbins,"Track Length (mm)","Counts")
    Hist1d("Y90Beta Track",Y90Distances,lenbins,"Track Length (mm)","Counts")

def PlotsforCs137(theEvents,title,lenbins=np.arange(1,80,10),enrg=500,ebins=np.arange(1,1000,60),trckthreshold=50,pltlimit=2,FileSave=False):
    Energys=[]
    cutE=[]
    Distances=[]
    FirstDistances=[]
    print(f"Number of Interacted Events is {len(theEvents)}")
    cnt=0
    for Event in theEvents:
        ElectronHits=theEvents[Event][0]
        Energys.append(theEvents[Event][1])

        for PID in ElectronHits:

            if(theEvents[Event][1]>enrg):
                title="Track with Cs137 " + str(theEvents[Event][1]) + " keV"
                Distances.append(ElectronHits[PID][1])
                FirstDistances.append(ElectronHits[PID][2])
                cutE.append(theEvents[Event][1])

                if(ElectronHits[PID][1]>trckthreshold):
                    if(cnt<=pltlimit):
                        Circle(ElectronHits[PID][0][1],ElectronHits[PID][0][2],title,FileSave)
                    cnt=cnt+1
    print (f" Number of Events after {enrg} keV energy cut is {len(Distances)}")
    print (f" Number of Events with {trckthreshold} mm and higher tracks is {cnt}")
    Hist1d("Cs137 Electron Energy Spectrum",Energys,ebins,"Energy (keV)","Counts")
    Hist1d("Cs137 Electron Energy Spectrum After Cut",cutE,ebins,"Energy (keV)","Counts")
    Hist1d("Cs137 First Hit Distance",FirstDistances,np.arange(1,200,10),"Length (mm)","Counts")
    Hist1d("Cs137 Electron ",Distances,lenbins,"Track Length (mm)","Counts")

def SaveFile(file,dic):
    np.save(file,dic,allow_pickle='TRUE')
def LoadNumpyFile(file):
    new_dict = np.load(file, allow_pickle='TRUE')
    print(new_dict.item())
    return new_dict


def main():

    # Get the Qualified Events
    #theEvents=RunEventsDic("/media/ilker/writable/DATA/Ba133_2mm_100000.h5",6) # For multiprocessing
    #theEvents=RunEventsDic("/media/ilker/writable/DATA/Cs137_2mm_100000.h5","/media/ilker/writable/DATA/Cs137_2mm",6) # For multiprocessing
    #theEvents=RunEventsDic("/media/ilker/writable/DATA/Cs137_2mm_1mm_away_1M_OCT5.h5","/media/ilker/writable/DATA/Cs137_2mm_1M",8,2000) # For multiprocessing
    theEvents=RunEventsDic("/media/ilker/writable/DATA/Ba133_1M.h5","/media/ilker/writable/DATA/Ba133_1M_v2",8,1000) # For multiprocessing
    #PlotEnergySpec(theEvents,"Cs137_2mm_1M_Oct")
    #PlotMotherEnergySpec(theEvents,"Gamma_Cs137_2mm_1M_Oct")


    #Plot RondomEvents From Region of Interest
    #PlotRandomEvents(theEvents,1,10,200,"Ba133_2mm")
    #PlotRandomEvents(theEvents,1,300,400,"Ba133_2mm")


    # Get the Qualified Events
    #theEvents=RunEventsDic("/media/ilker/writable/DATA/Cs137_10mm_100000.h5",6) # For multiprocessing
    #theEvents=RunEventsDic("/home/ilker/Dropbox/nexus/build2/source/Cs137_2mm1M.h5",8) # For multiprocessing
    #PlotEnergySpec(theEvents,"Energy Spectrum of Cs137_10mm_100k")

    #Plot RondomEvents From Region of Interest
    #PlotRandomEvents(theEvents,2,10,200,"Cs137_10mm_100k")
    #PlotRandomEvents(theEvents,2,300,400,"Cs137_10mm_100k")




if __name__ == "__main__":
    st=time.time()
    main()
    print("--- Finished in  %s minutes ---" %((time.time()-st)/60))
