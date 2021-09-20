import numpy             as np
import matplotlib.pyplot as plt
import os
import glob
import h5py as h
import multiprocessing as mtp
import time
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

    #plt.savefig('/Users/austinmcdonald/Desktop/'+title + '.png')
    plt.show()
def Plot2d(title,x,y):
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
def getEnergys(file):
    data = h.File(file,'r')
    len1=data['MC']['configuration'][2][-1]
    Energys = []
    for x in range(0,int(len1)):
        Current_Event = x
        Current_Hit_Mask = data["MC"]['hits']['event_id'] == Current_Event

        Current_Particle_Mask = data["MC"]['particles']['event_id'] == Current_Event

        Electron_Mask = data["MC"]['particles'][Current_Particle_Mask]['particle_name'] == b'e-'

        Electron_PIDS = data["MC"]['particles'][Current_Particle_Mask][Electron_Mask]['particle_id']

        A = data["MC"]['hits'][Current_Hit_Mask]['particle_id']
        Hit_Electron_Maks = np.in1d(A, Electron_PIDS)

        MYE = data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]['energy'].sum()*1e3
        if (MYE>450):
            print("Energy {} \t event {} ".format(MYE,Current_Event))
        Energys.append(data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]['energy'].sum()*1e3)
    Energys = np.array(Energys)
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

    '''for i in ElectronData["particle_id"]:
        PID_mask=ElectronData["particle_id"]==i
        count = np.count_nonzero(PID_mask)
        PIDWithCounts[count]=PID_mask

    max_key=max(list(PIDWithCounts.keys()))
    ElectronsWithMaxTrack=ElectronData[PIDWithCounts[max_key]]
    '''
    Xhits=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]["x"]
    Yhits=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]["y"]
    Zhits=data["MC"]['hits'][Current_Hit_Mask][Hit_Electron_Maks]["z"]
    Tracks=np.array([Xhits,Yhits,Zhits])

    return Tracks


def fun(val):
    print(f'hello {val}')
    time.sleep(val)
    print (f'Finishing time with val {val}')
    print('parent process id:',os.getppid())
    print('process id:',os.getpid())

def main():

    p1 = mtp.Process(target=fun, args=(3,))
    p1.start()
    p1.join()
    print(f'Process p1 is alive: {p1.is_alive()}')
    p2 = mtp.Process(target=fun, args=(3,))
    p2.start()
    p2.join()
    p3 = mtp.Process(target=fun, args=(3,))
    p3.start()
    p3.join()

if __name__ == "__main__":
    #p=mtp.Process(target=getEnergys,args=("/home/ilker/Dropbox/nexus/build/source/Cs137_2mm12.h5"))
    print('starting Main')
    main()
    print('starting Main')
    main()
    #Hist1d('Energy Spectrum of Cs137',p,np.arange(1,800,10),0,0)