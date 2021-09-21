import multiprocessing as mtp
import numpy as np
import matplotlib.pyplot as plt
import h5py as h
import time
import os
import random
#Create a custom 2d Histogram
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

    plt.savefig('/home/ilker/Dropbox/nexus/build/source/'+title + '.png')
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
        savetit='/home/ilker/Pictures/Sim/'+title + '.png'
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
        savetit='/home/ilker/Pictures/Sim/'+title + '.png'
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
def RunEventsDic(file,Processes=4,TotalLimit=0):
    data=ReadFile(file)
    theEvents={}
    if TotalLimit==0:
        TotalEvents=int(data['MC']['configuration'][2][1])
    else:
        TotalEvents=TotalLimit
    print(f"Processing Your Events at {file}")
    print(f"Currently there are {TotalEvents} events")

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

    return theEvents



def EnergyAndTracks(data,xRLower,xRHigher,theEvents,EventLimit=0,fudicalR=35,fdcount=True,AllEvents=False):
    TotalEvents=xRHigher
    if(EventLimit>0):
        TotalEvents=EventLimit

    FudicalCount=0
    for Current_Event in range(xRLower,TotalEvents):

        Current_Hit_Mask = data['MC']['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event

        Electron_Mask = data['MC']['particles'][Current_Particle_Mask]['particle_name'] == b'e-'

        #Gamma_Mask = data['MC']['particles'][Current_Particle_Mask]['particle_name'] == b'gamma'
        #Gamma_PIDS = data['MC']['particles'][Current_Particle_Mask]['particle_id']

        #Electron_MIDS = data["MC"]['particles'][Current_Particle_Mask][Electron_Mask]['mother_id']

        #Mothers=np.in1d(Gamma_PIDS,Electron_MIDS)

        #MotherEnergy=data['MC']['particles'][Current_Particle_Mask][Mothers]['kin_energy']
        #print(MotherEnergy*1000)

        Electron_PIDS = data["MC"]['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']

        A = data["MC"]['hits'][Current_Hit_Mask]['particle_id']
        Hit_Electron_Maks = np.in1d(A, Electron_PIDS)
        ElectronData=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]

        TotalEventEnergy=ElectronData["energy"].sum()*1e3
        Xhits=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]["x"]
        Yhits=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]["y"]
        Zhits=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]["z"]

        #Checking the Radius if in the boundaries
        fRadius=np.sqrt(Yhits*Yhits+Zhits*Zhits)
        fRadius_Mask=fRadius<=fudicalR

        TrksCntIntheFudi=np.count_nonzero(fRadius_Mask)
        if(not AllEvents):
            if((TrksCntIntheFudi==0)):
                continue
        FudicalCount=FudicalCount+1

        Tracks=np.array([Xhits,Yhits,Zhits])
        theEvents[Current_Event]=[Tracks,TotalEventEnergy,TrksCntIntheFudi]


    if(fdcount):
        print(f" There are only {FudicalCount} out of {TotalEvents} in the fudical volume")


def PlotRandomEvents(theEvents,NPlots,LowE,HighE,pretitle,FileSave=True):

    QualfEnergys=[]
    for Event in theEvents:
        if(theEvents[Event][1]>=LowE and theEvents[Event][1]<=HighE):
            QualfEnergys.append(Event)
    #print(len(QualfEnergys))
    for plot in range(0,NPlots):
        Event=random.choice(QualfEnergys)
        title=str(pretitle) + "_Tracks_" + str(round(theEvents[Event][1],2)) + "keV"
        Circle(theEvents[Event][0][1],theEvents[Event][0][2],title,FileSave)


def PlotEnergySpec(theEvents,title,binss=np.arange(1,500,10),FileSave=True):
    Energys=[]
    for Event in theEvents:
        Energys.append(theEvents[Event][1])
    Hist1d(title,Energys,binss,0,0,FileSave)
def main():

    # Get the Qualified Events
    theEvents=RunEventsDic("/home/ilker/Dropbox/nexus/build2/source/Ba133_1M.h5",8) # For multiprocessing
    #theEvents=RunEventsDic("/home/ilker/Dropbox/nexus/build2/source/Cs137_2mm1M.h5",8) # For multiprocessing
    PlotEnergySpec(theEvents,"Energy Spectrum of Ba133_2mm")

    #Plot RondomEvents From Region of Interest
    PlotRandomEvents(theEvents,5,100,200,"Ba133")
    PlotRandomEvents(theEvents,5,400,500,"Ba133")

if __name__ == "__main__":
    st=time.time()
    main()
    print("--- Finished in  %s minutes ---" %((time.time()-st)/60))
