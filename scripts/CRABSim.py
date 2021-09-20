import multiprocessing as mtp
import numpy as np
import matplotlib.pyplot as plt
import h5py as h
import time
import os
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

def Hist1d(title,Energys,Bins,xlimit,ylimit,limits=False):
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

    plt.savefig('/home/ilker/Pictures/Sim/'+title + '.png')


def ReadFile(file):
    data = h.File(file,'r')
    return data
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
def Circle(pltX,pltY,title):
    fig,ax = plt.subplots(figsize=(8,8))
    plt.plot(pltX,pltY)
    plt.scatter(pltX,pltY)
    plt.title(title)
    circle1=plt.Circle((70,70),0.1,color='r',fill=False)
    ax.add_patch(circle1)
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.ylabel("Z mm",fontsize=26)
    plt.xlabel("Y mm",fontsize=26)
    plt.tick_params('both', length=10, width=2, which='major')
    plt.tick_params('both', length=5, width=1, which='minor')
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def RunEvents(file,Processes=4,AnaFast=True):
    data=ReadFile(file)
    Energys=[]
    TotalEvents=int(data['MC']['configuration'][2][1])
    print(f"Currently there are {TotalEvents}")
    print(f"Processing Your Events at {file}")
    if(AnaFast):
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

def main():

    #EnergyValues=RunEvents("/home/ilker/Dropbox/nexus/build/source/Ba133_2mm222.h5") # For multiprocessing
    EnergyValues=RunEvents("/home/ilker/Dropbox/nexus/build2/source/Cs137_2mm1M.h5",6) # For multiprocessing
    #EnergyValues=RunEvents("/home/ilker/Dropbox/nexus/build/source/Ba133_2mm222.h5",False)

    Hist1d("Cs137_2mm_1M",EnergyValues,np.arange(1,500,10),0,0)
if __name__ == "__main__":
    st=time.time()
    main()
    print("--- %s minutes ---" %((time.time()-st)/60))
