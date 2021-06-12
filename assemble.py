# assemble.py
import logging

def energy(Jobs,level):
    #loMonJobs=Jobs[0]
    #loPairJobs=Jobs[1]
    #loTrimerJobs=Jobs[2]
    #clusterJob=Jobs[3]
    #hiMonJobs=Jobs[4]
    #hiPairJobs=Jobs[5]
    #hiTrimerJobs=Jobs[6]
    clusterJob=Jobs[0]
    if level>1:
        loMonJobs=Jobs[1]
        hiMonJobs=Jobs[2]
        loPairJobs=Jobs[3]
        hiPairJobs=Jobs[4]
        if level>2:
            loTrimerJobs=Jobs[5]
            hiTrimerJobs=Jobs[6]


    E=clusterJob.energy
    try:
        n=len(loMonJobs)
    except:
        n=1

    # Exit if only a cluster comp. 
    if level==1:
        return E
    if level==2:
        oneB_factor=n-2
        twoB_factor=-1
    if level==3:
        oneB_factor=(n-2)*(n-3)/-2.0
        twoB_factor=n-3

    for i in loMonJobs:
        #E -= (n-2)*(n-3)/2.0*i.energy
        E += oneB_factor*i.energy
    for i in hiMonJobs:
        #E += (n-2)*(n-3)/2.0*i.energy
        E -= oneB_factor*i.energy

    for i in loPairJobs:
        #E += (n-3)*i.energy
        E += twoB_factor*i.energy
    for i in hiPairJobs:
       # E -= (n-3)*i.energy
        E -= twoB_factor*i.energy

    # Exit if 2-body
    if level==2:
        return E

    for i in loTrimerJobs:
        E -= i.energy
    for i in hiTrimerJobs:
        E += i.energy

    return E


def gradient(Jobs,level):
    #loMonJobs=Jobs[0]
    #loPairJobs=Jobs[1]
    #loTrimerJobs=Jobs[2]
    #clusterJob=Jobs[3]
    #hiMonJobs=Jobs[4]
    #hiPairJobs=Jobs[5]
    #hiTrimerJobs=Jobs[6]
    clusterJob=Jobs[0]
    if level>1:
        loMonJobs=Jobs[1]
        hiMonJobs=Jobs[2]
        loPairJobs=Jobs[3]
        hiPairJobs=Jobs[4]

        if level>2:
            loTrimerJobs=Jobs[5]
            hiTrimerJobs=Jobs[6]

    G = clusterJob.rgrad
    natom=clusterJob.natom
    try:
        n=len(loMonJobs)
    except:
        n=1

    if level==1:
        return G

    if level==2:
        oneB_factor=n-2
        twoB_factor=-1

    if level==3:
        oneB_factor=(n-2)*(n-3)/-2.0
        print("oneB factor = ",oneB_factor)
        twoB_factor=n-3


    for i in loMonJobs:
        a=0
        for j in i.atoms:
            if j==0:
                print("Before Subtracting from atom 0")
                print(G[0])
            for k in range(3):
                #G[j][k]-=(n-2)*(n-3)/2.0*i.rgrad[a][k]
                G[j][k]+=oneB_factor*i.rgrad[a][k]
            if j==0:
                print("After subtraction atom 0 = ",G[0])
            a+=1
            
    for i in hiMonJobs:
        a=0
        for j in i.atoms:
            if j==0:
                print("Adding to atom 0")
            for k in range(3):
                #G[j][k]+=(n-2)*(n-3)/2.0*i.rgrad[a][k]
                G[j][k]-=oneB_factor*i.rgrad[a][k]
            if j==0:
                print("After addition atom 0 = ",G[0])
            a+=1
    print("Coleman")
    print("Update after mons")
    print(G)
     
    for i in loPairJobs:
        a=0
        for j in i.atoms:
            for k in range(3):
                #G[j][k]+=(n-3)*i.rgrad[a][k]
                G[j][k]+=twoB_factor*i.rgrad[a][k]
            a+=1

    for i in hiPairJobs:
        a=0
        for j in i.atoms:
            for k in range(3):
                #G[j][k]-=(n-3)*i.rgrad[a][k]
                G[j][k]-=twoB_factor*i.rgrad[a][k]
            a+=1
    print("Coleman")
    print("Update after pairs")
    print(G)
   
    # Exit if 2-body
    if level==2:
        return G

    for i in loTrimerJobs:
        a=0
        for j in i.atoms:
            for k in range(3):
                G[j][k]-=i.rgrad[a][k]
            a+=1

    for i in hiTrimerJobs:
        a=0
        for j in i.atoms:
            for k in range(3):
                G[j][k]+=i.rgrad[a][k]
            a+=1

    print("Coleman")
    print("Update after trimers")
    print(G)


    return G

def hessian(Jobs,level):
    asslogger=logging.getLogger('Hessian Assembly')
    #loMonJobs=Jobs[0]
    #loPairJobs=Jobs[1]
    #loTrimerJobs=Jobs[2]
    #clusterJob=Jobs[3]
    #hiMonJobs=Jobs[4]
    #hiPairJobs=Jobs[5]
    #hiTrimerJobs=Jobs[6]
    clusterJob=Jobs[0]
    if level>1:
        loMonJobs=Jobs[1]
        hiMonJobs=Jobs[2]
        loPairJobs=Jobs[3]
        hiPairJobs=Jobs[4]

        if level>2:
            loTrimerJobs=Jobs[5]
            hiTrimerJobs=Jobs[6]

    H=clusterJob.hessian

    natom=clusterJob.natom

    try:
        n=len(loMonJobs)
    except:
        n=1

    if level==1:
        return H
    if level==2:
        oneB_factor=n-2
        twoB_factor=-1
    if level==3:
        oneB_factor=(n-2)*(n-3)/-2.0
        twoB_factor=n-3


    for mons in loMonJobs:
        a=0
        for i in mons.atoms:
            for j in range(3):
                b=0
                for k in mons.atoms:
                    for l in range(3):
                      #  H[i][j][k][l]-=(n-2)*(n-3)/2.0*mons.hessian[a][j][b][l]
                        H[i][j][k][l]+=oneB_factor*mons.hessian[a][j][b][l]
                    b+=1
            a+=1

    for mons in hiMonJobs:
        a=0
        for i in mons.atoms:
            for j in range(3):
                b=0
                for k in mons.atoms:
                    for l in range(3):
                        #H[i][j][k][l]+=(n-2)*(n-3)/2.0*mons.hessian[a][j][b][l]
                        H[i][j][k][l]-=oneB_factor*mons.hessian[a][j][b][l]
                    b+=1
            a+=1


    for pairs in loPairJobs:
        a=0
        for i in pairs.atoms:
            for j in range(3):
                b=0
                for k in pairs.atoms:
                    for l in range(3):
                       # H[i][j][k][l]+=(n-3)*pairs.hessian[a][j][b][l]
                        H[i][j][k][l]+=twoB_factor*pairs.hessian[a][j][b][l]
                    b+=1
            a+=1


    for pairs in hiPairJobs:
        a=0
        for i in pairs.atoms:
            for j in range(3):
                b=0
                for k in pairs.atoms:
                    for l in range(3):
                        #H[i][j][k][l]-=(n-3)*pairs.hessian[a][j][b][l]
                        H[i][j][k][l]-=twoB_factor*pairs.hessian[a][j][b][l]
                    b+=1
            a+=1

    if level==2:
        return H

    for trimers in loTrimerJobs:
        a=0
        for i in trimers.atoms:
            for j in range(3):
                b=0
                for k in trimers.atoms:
                    for l in range(3):
                        H[i][j][k][l]-=trimers.hessian[a][j][b][l]
                    b+=1
            a+=1

    for trimers in hiTrimerJobs:
        a=0
        for i in trimers.atoms:
            for j in range(3):
                b=0
                for k in trimers.atoms:
                    for l in range(3):
                        H[i][j][k][l]+=trimers.hessian[a][j][b][l]
                    b+=1
            a+=1

    asslogger.info('\nHessian after all terms\n%s' %H)
    return H


