#newass.py
import logging,sys
import numpy as np

def newenergy(Jobs,level):
    al=logging.getLogger('Energy Assembly')
    clusterJob=Jobs[0]
    if level==0:
        return clusterJob.energy
    if level>0:
        loMonJobs=Jobs[1]
        hiMonJobs=Jobs[2]

        if level>1:
            loPairJobs=Jobs[3]
            hiPairJobs=Jobs[4]

            if level>2:
                loTrimerJobs=Jobs[5]
                hiTrimerJobs=Jobs[6]

                if level>3:
                    loTetradJobs=Jobs[7]
                    hiTetradJobs=Jobs[8]

    # initialize counters for debugging
    monCount=0
    pairCount=0
    triCount=0
    tetCount=0
    Eref=clusterJob.energy

    E=0.0

    if level>0:
    # Add E_1 from hi-level
        for i in hiMonJobs:
            #al.info('Checking this mon:%s',i.atoms)
            E+=i.energy # Add this E
            monCount+=1
        al.info('Calculated E1 = %s',E)
        error=(E-Eref)*627.51
        al.info('Error = %s',error)

    # Add E_2 from hi-level
    if level>1:
        for i in hiPairJobs:
            E+=i.energy
            pairCount+=1 # Add dimer energies
            # Subtract monomer overlap
            for j in hiMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    #al.info("Subtracting out this mon: %s",j.atoms)
                    E-=j.energy

        al.info('Calculated E2 = %s', E)
        error=(E-Eref)*627.51
        al.info('Error = %s',error)

    if level>2:
        # Add E_3 from hi-level
        for i in hiTrimerJobs:
            E+=i.energy
            triCount+=1
            # Subtract pair overlap
            for j in hiPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E-=j.energy
                    # If match, add back in mon overlap
                    for k in hiMonJobs:
                        monMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            E+=k.energy
            # Subtract mon overlap
            for j in hiMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E-=j.energy
        al.info("Calculated E3 = %s",E)
        error=(E-Eref)*627.51
        al.info('Error = %s',error)

    # Add E4 from hi level
    if level>3:
        for i in hiTetradJobs:
            E+=i.energy
            # Subtract trimer overlap
            for j in hiTrimerJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E-=j.energy
                    # Add back in pair overlap
                    for k in hiPairJobs:
                        pMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                pMatch=False
                                break
                        if pMatch:
                            E+=k.energy
                            # Subtract out mon overlap from pair
                            for l in hiMonJobs:
                                mMatch=True
                                for atom in l.atoms:
                                    if atom in k.atoms:
                                        pass
                                    else:
                                        mMatch=False
                                        break
                                if mMatch:
                                    E-=l.energy

                    # Add back in mon overlap from trimer 
                    for k in hiMonJobs:
                        mMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                mMatch=False
                                break
                        if mMatch:
                            E+=k.energy

            # Subtract pair overlap from tetrad 
            for j in hiPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E-=j.energy
                    # Add back in mon overlap from pair
                    for k in hiMonJobs:
                        monMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            E+=k.energy

            # Subtract mon overlap from tetrad
            for j in hiMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E-=j.energy


        al.info("Calculated E4 = %s",E)
        error=(E-Eref)*627.51
        al.info('Error = %s',error)

    # Now copy and paste hi-level shit and change to lo-level

    # Subtract low level E1
    if level>0:
        for i in loMonJobs:
            #al.info('Checking this mon:%s',i.atoms)
            E-=i.energy # Add this E
            monCount+=1

    # Subtract low level E2
    if level>1:
        for i in loPairJobs:
            E-=i.energy
            pairCount+=1 # Add dimer energies
            # Add monomer overlap
            for j in loMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    #al.info("Subtracting out this mon: %s",j.atoms)
                    E+=j.energy


    if level>2:
        # Subtract low level E3
        for i in loTrimerJobs:
            E-=i.energy
            triCount+=1
            # Add pair overlap
            for j in loPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E+=j.energy
                    # If match, sub out mon overlap
                    for k in loMonJobs:
                        monMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            E-=k.energy
            # Add mon overlap
            for j in loMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E+=j.energy

    # Subtract low E4 
    if level>3:
        for i in loTetradJobs:
            E-=i.energy
            # Add trimer overlap
            for j in loTrimerJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E+=j.energy
                    # Subtract out pair overlap
                    for k in loPairJobs:
                        pMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                pMatch=False
                                break
                        if pMatch:
                            E-=k.energy
                            # Add mon overlap from pair
                            for l in loMonJobs:
                                mMatch=True
                                for atom in l.atoms:
                                    if atom in k.atoms:
                                        pass
                                    else:
                                        mMatch=False
                                        break
                                if mMatch:
                                    E+=l.energy

                    # Sub out mon overlap from trimer 
                    for k in loMonJobs:
                        mMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                mMatch=False
                                break
                        if mMatch:
                            E-=k.energy

            # Add pair overlap from tetrad 
            for j in loPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E+=j.energy
                    # Sub out mon overlap from pair
                    for k in loMonJobs:
                        monMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            E-=k.energy

            # Add mon overlap from tetrad
            for j in loMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    E+=j.energy



    # Now Add cluster Energy to get n-body:many-body energy
    E+=Eref

    al.info('Calculated E = %s',E)
    error=(E-Eref)*627.51
    al.info('Difference from low level (kcal/mol) = %s',error)

    return E

def newgrad(Jobs,level):
    al=logging.getLogger('Gradient Assembly')
    clusterJob=Jobs[0]
    if level==0:
        return clusterJob.rgrad
    if level>0:
        loMonJobs=Jobs[1]
        hiMonJobs=Jobs[2]

        if level>1:
            loPairJobs=Jobs[3]
            hiPairJobs=Jobs[4]

            if level>2:
                loTrimerJobs=Jobs[5]
                hiTrimerJobs=Jobs[6]

                if level>3:
                    loTetradJobs=Jobs[7]
                    hiTetradJobs=Jobs[8]

    G=clusterJob.rgrad

    if level>0:
    # Add G_1 from hi-level
        for i in hiMonJobs:
            a=0
            for j in i.atoms:
                for k in range(3):
                    G[j][k]+=i.rgrad[a][k]
                a+=1

    # Add G_2 from hi-level
    if level>1:
        for i in hiPairJobs:
            a=0
            for j in i.atoms:
                for k in range(3):
                    G[j][k]+=i.rgrad[a][k]
                a+=1

            # Subtract monomer overlap from hiPair_i
            for j in hiMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    # Hit - subtract the grad
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]-=j.rgrad[a][k]
                        a+=1
    # Add G_3 from hi-level
    if level>2:
        for i in hiTrimerJobs:
            a=0
            for j in i.atoms:
                for k in range(3):
                    G[j][k]+=i.rgrad[a][k]
                a+=1

            # Subtract Pair overlap
            for j in hiPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    # Hit - subtract the pair grad
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]-=j.rgrad[a][k]
                        a+=1
                    # Must add back mon grad if matched
                    for k in hiMonJobs:
                        monMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            # Hit - add back mon from pair overlap
                            for atom in k.atoms:
                                for l in range(3):
                                    G[atom][l]+=k.rgrad[a][l]
                                a+=1

            # Subtract mon overlap (from trimer)
            for j in hiMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    # Hit - subtract the mon grad
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]-=j.rgrad[a][k]
                        a+=1


    # Add G4 from hi-level
    if level>3:
        for i in hiTetradJobs:
            a=0
            for atom in i.atoms:
                for k in range(3):
                    G[atom][k]+=i.rgrad[a][k]
                a+=1

            # Subtract trimer overlap from tetrad
            for j in hiTrimerJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]-=j.rgrad[a][k]
                        a+=1

                    # Add back pair overlap from trimer
                    for l in hiPairJobs:
                        pMatch=True
                        for atom in l.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                pMatch=False
                                break
                        if pMatch:
                            a=0
                            for atom in l.atoms:
                                for k in range(3):
                                    G[atom][k]+=l.rgrad[a][k]
                                a+=1

                            # Sub back out mon from pair
                            for m in hiMonJobs:
                                mMatch=True
                                for atom in m.atoms:
                                    if atom in l.atoms:
                                        pass
                                    else:
                                        mMatch=False
                                        break
                                if mMatch:
                                    a=0
                                    for atom in m.atoms:
                                        for k in range(3):
                                            G[atom][k]-=m.rgrad[a][k]
                                        a+=1

                    # Add back mon overlap from trimer
                    for l in hiMonJobs:
                        mMatch=True
                        for atom in l.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                mMatch=False
                                break
                        if mMatch:
                            a=0
                            for atom in l.atoms:
                                for k in range(3):
                                    G[atom][k]+=l.rgrad[a][k]
                                a+=1
            # Subtract pair overlap from tetrad
            for j in hiPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]-=j.rgrad[a][k]
                        a+=1

                    # Add back mon overlap from pair
                    for l in hiMonJobs:
                        mMatch=True
                        for atom in l.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                mMatch=False
                                break
                        if mMatch:
                            a=0
                            for atom in l.atoms:
                                for k in range(3):
                                    G[atom][k]+=l.rgrad[a][k]
                                a+=1

            # Subtract mon overlap from tetrad
            for j in hiMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]-=j.rgrad[a][k]
                        a+=1
##############################################
#       Subtracting off low level terms
###############################################

    # Subtract G_1
    if level>0:
        for i in loMonJobs:
            a=0
            for atom in i.atoms:
                for k in range(3):
                    G[atom][k]-=i.rgrad[a][k]
                a+=1

    # Subtract G_2
    if level>1:
        for i in loPairJobs:
            a=0
            for atom in i.atoms:
                for k in range(3):
                    G[atom][k]-=i.rgrad[a][k]
                a+=1

            # Add back in mon overlap from loPair
            for j in loMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]+=j.rgrad[a][k]
                        a+=1

    # Subtract G_3
    if level>2:
        for i in loTrimerJobs:
            a=0
            for atom in i.atoms:
                for k in range(3):
                    G[atom][k]-=i.rgrad[a][k]
                a+=1

            # Add back pair overlap from loTrimer
            for j in loPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]+=j.rgrad[a][k]
                        a+=1

                    # Must subtract out mon overlap of pair
                    for k in loMonJobs:
                        monMatch=True
                        for atom in k.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for atom in k.atoms:
                                for l in range(3):
                                    G[atom][l]-=k.rgrad[a][l]
                                a+=1

            # Add back mon overlap from loTrimer
            for j in loMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]+=j.rgrad[a][k]
                        a+=1


    # Sub G4 from lo-level
    if level>3:
        for i in loTetradJobs:
            a=0
            for atom in i.atoms:
                for k in range(3):
                    G[atom][k]-=i.rgrad[a][k]
                a+=1

            # Add trimer overlap from tetrad
            for j in loTrimerJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]+=j.rgrad[a][k]
                        a+=1

                    # Sub back pair overlap from trimer
                    for l in loPairJobs:
                        pMatch=True
                        for atom in l.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                pMatch=False
                                break
                        if pMatch:
                            a=0
                            for atom in l.atoms:
                                for k in range(3):
                                    G[atom][k]-=l.rgrad[a][k]
                                a+=1

                            # Add back in mon from pair
                            for m in loMonJobs:
                                mMatch=True
                                for atom in m.atoms:
                                    if atom in l.atoms:
                                        pass
                                    else:
                                        mMatch=False
                                        break
                                if mMatch:
                                    a=0
                                    for atom in m.atoms:
                                        for k in range(3):
                                            G[atom][k]+=m.rgrad[a][k]
                                        a+=1

                    # Sub back mon overlap from trimer
                    for l in loMonJobs:
                        mMatch=True
                        for atom in l.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                mMatch=False
                                break
                        if mMatch:
                            a=0
                            for atom in l.atoms:
                                for k in range(3):
                                    G[atom][k]-=l.rgrad[a][k]
                                a+=1
            # Add pair overlap from tetrad
            for j in loPairJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]+=j.rgrad[a][k]
                        a+=1

                    # Sub back mon overlap from pair
                    for l in loMonJobs:
                        mMatch=True
                        for atom in l.atoms:
                            if atom in j.atoms:
                                pass
                            else:
                                mMatch=False
                                break
                        if mMatch:
                            a=0
                            for atom in l.atoms:
                                for k in range(3):
                                    G[atom][k]-=l.rgrad[a][k]
                                a+=1

            # Add mon overlap from tetrad
            for j in loMonJobs:
                match=True
                for atom in j.atoms:
                    if atom in i.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for atom in j.atoms:
                        for k in range(3):
                            G[atom][k]+=j.rgrad[a][k]
                        a+=1
    return G

def newhess(Jobs,level):
    al=logging.getLogger('Hessian Assembly')
    clusterJob=Jobs[0]
    if level==0:
        return clusterJob.hessian
    if level>0:
        loMonJobs=Jobs[1]
        hiMonJobs=Jobs[2]

        if level>1:
            loPairJobs=Jobs[3]
            hiPairJobs=Jobs[4]

            if level>2:
                loTrimerJobs=Jobs[5]
                hiTrimerJobs=Jobs[6]

                if level>3:
                    loTetradJobs=Jobs[7]
                    hiTetradJobs=Jobs[8]

    H=clusterJob.hessian

    # Add H_1 from hi-level
    if level>0: # Line not necessary, pure aesthetics
        for mon in hiMonJobs:
            a=0
            for i in mon.atoms:
                for j in range(3):
                    b=0
                    for k in mon.atoms:
                        for l in range(3):
                            H[i][j][k][l]+=mon.hessian[a][j][b][l]
                        b+=1
                a+=1

    # Add H_2 from hi-level
    if level>1:
        for pair in hiPairJobs:
            a=0
            for i in pair.atoms:
                for j in range(3):
                    b=0
                    for k in pair.atoms:
                        for l in range(3):
                            H[i][j][k][l]+=pair.hessian[a][j][b][l]
                        b+=1
                a+=1



            # Sub hi-mon overlap from pair
            for mon in hiMonJobs:
                match=True
                for atom in mon.atoms:
                    if atom in pair.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in mon.atoms:
                        for j in range(3):
                            b=0
                            for k in mon.atoms:
                                for l in range(3):
                                    H[i][j][k][l]\
                                    -=mon.hessian[a][j][b][l]
                                b+=1
                        a+=1

    # Add H_3 from hi-level
    if level>2:
        for trimer in hiTrimerJobs:
            a=0
            for i in trimer.atoms:
                for j in range(3):
                    b=0
                    for k in trimer.atoms:
                        for l in range(3):
                            H[i][j][k][l]+=trimer.hessian[a][j][b][l]
                        b+=1
                a+=1

            # Sub pair overlap from trimer
            for pair in hiPairJobs:
                match=True
                for atom in pair.atoms:
                    if atom in trimer.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in pair.atoms:
                        for j in range(3):
                            b=0
                            for k in pair.atoms:
                                for l in range(3):
                                    H[i][j][k][l]-=\
                                        pair.hessian[a][j][b][l]
                                b+=1
                        a+=1
                        
                # Add back mon overlap from pair
                    for mon in hiMonJobs:
                        monMatch=True
                        for atom in mon.atoms:
                            if atom in pair.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for i in mon.atoms:
                                for j in range(3):
                                    b=0
                                    for k in mon.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]+=\
                                             mon.hessian[a][j][b][l]
                                        b+=1
                                a+=1

            # Subtract mon overlap from trimer
            for mon in hiMonJobs:
                match=True
                for atom in mon.atoms:
                    if atom in trimer.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in mon.atoms:
                        for j in range(3):
                            b=0
                            for k in mon.atoms:
                                for l in range(3):
                                    H[i][j][k][l]-=\
                                     mon.hessian[a][j][b][l]

                                b+=1
                        a+=1

    # Add H_4 from hi-level
    if level>3:
        for tetrad in hiTetradJobs:
            a=0
            for i in tetrad.atoms:
                for j in range(3):
                    b=0
                    for k in tetrad.atoms:
                        for l in range(3):
                            H[i][j][k][l]+=\
                                tetrad.hessian[a][j][b][l]

                        b+=1
                a+=1

            # Sub trimer overlap from tetrad
            for trimer in hiTrimerJobs:
                match=True
                for atom in trimer.atoms:
                    if atom in tetrad.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in trimer.atoms:
                        for j in range(3):
                            b=0
                            for k in trimer.atoms:
                                for l in range(3):
                                    H[i][j][k][l]-=\
                                      trimer.hessian[a][j][b][l]
                                b+=1
                        a+=1

                    # Add pair overlap from trimer
                    for pair in hiPairJobs:
                        pairMatch=True
                        for atom in pair.atoms:
                            if atom in trimer.atoms:
                                pass
                            else:
                                pairMatch=False
                                break
                        if pairMatch:
                            a=0
                            for i in pair.atoms:
                                for j in range(3):
                                    b=0
                                    for k in pair.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]+=\
                                             pair.hessian[a][j][b][l]

                                        b+=1
                                a+=1

                            # Sub mon overlap from pair
                            for mon in hiMonJobs:
                                monMatch=True
                                for atom in mon.atoms:
                                    if atom in pair.atoms:
                                        pass
                                    else:
                                        monMatch=False
                                        break
                                if monMatch:
                                    a=0
                                    for i in mon.atoms:
                                        for j in range(3):
                                            b=0
                                            for k in mon.atoms:
                                                for l in range(3):
                                                    H[i][j][k][l]-=\
                                                 mon.hessian[a][j][b][l]
                                                b+=1
                                        a+=1

                    # Add mon overlap from trimer
                    for mon in hiMonJobs:
                        monMatch=True
                        for atom in mon.atoms:
                            if atom in trimer.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for i in mon.atoms:
                                for j in range(3):
                                    b=0
                                    for k in mon.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]+=\
                                             mon.hessian[a][j][b][l]
                                        b+=1
                                a+=1

            # Sub pair overlap from tetrad
            for pair in hiPairJobs:
                match=True
                for atom in pair.atoms:
                    if atom in tetrad.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in pair.atoms:
                        for j in range(3):
                            b=0
                            for k in pair.atoms:
                                for l in range(3):
                                    H[i][j][k][l]-=\
                                        pair.hessian[a][j][b][l]
                                b+=1
                        a+=1

                    # Add mon overlap from pair
                    for mon in hiMonJobs:
                        monMatch=True
                        for atom in mon.atoms:
                            if atom in pair.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for i in mon.atoms:
                                for j in range(3):
                                    b=0
                                    for k in mon.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]+=\
                                             mon.hessian[a][j][b][l]
                                        b+=1
                                a+=1

                                    
            # Sub mon overlap from tetrad
            for mon in hiMonJobs:
                match=True
                for atom in mon.atoms:
                    if atom in tetrad.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in mon.atoms:
                        for j in range(3):
                            b=0
                            for k in mon.atoms:
                                for l in range(3):
                                    H[i][j][k][l]-=\
                                      mon.hessian[a][j][b][l]
                                b+=1
                        a+=1

#############################################################
#     Low Level Code        
#
#############################################################

    # Subtract H_1(low) and H_2(low)
    for mon in loMonJobs:
        a=0
        for i in mon.atoms:
            for j in range(3):
                b=0
                for k in mon.atoms:
                    for l in range(3):
                        H[i][j][k][l]-=mon.hessian[a][j][b][l]
                    b+=1
            a+=1

    # Sub H_2
    if level>1:
        for pair in loPairJobs:
            a=0
            for i in pair.atoms:
                for j in range(3):
                    b=0
                    for k in pair.atoms:
                        for l in range(3):
                            H[i][j][k][l]-=pair.hessian[a][j][b][l]
                        b+=1
                a+=1

            # Add back loMon overlap from pair
            for mon in loMonJobs:
                match=True
                for atom in mon.atoms:
                    if atom in pair.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in mon.atoms:
                        for j in range(3):
                            b=0
                            for k in mon.atoms:
                                for l in range(3):
                                    H[i][j][k][l]+=\
                                        mon.hessian[a][j][b][l]
                                b+=1
                        a+=1

    # Sub H_3
    if level>2:
        for trimer in loTrimerJobs:
            a=0
            for i in trimer.atoms:
                for j in range(3):
                    b=0
                    for k in trimer.atoms:
                        for l in range(3):
                            H[i][j][k][l]-=trimer.hessian[a][j][b][l]
                        b+=1
                a+=1

        # Add loPair overlap from trimer
            for pair in loPairJobs:
                match=True
                for i in pair.atoms:
                    if i in trimer.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in pair.atoms:
                        for j in range(3):
                            b=0
                            for k in pair.atoms:
                                for l in range(3):
                                    H[i][j][k][l]+=\
                                      pair.hessian[a][j][b][l]
                                b+=1
                        a+=1

                    # Sub mon overlap from pair    
                    for mon in loMonJobs:
                        monMatch=True
                        for atom in mon.atoms:
                            if atom in pair.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for i in mon.atoms:
                                for j in range(3):
                                    b=0
                                    for k in mon.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]-=\
                                             mon.hessian[a][j][b][l]
                                        b+=1
                                a+=1

            # Add mon overlap from trimer
            for mon in loMonJobs:
                match=True
                for atom in mon.atoms:
                    if atom in trimer.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in mon.atoms:
                        for j in range(3):
                            b=0
                            for k in mon.atoms:
                                for l in range(3):
                                    H[i][j][k][l]+=\
                                        mon.hessian[a][j][b][l]
                                b+=1
                        a+=1

    # Subtract H_4 from lo-level
    if level>3:
        for tetrad in loTetradJobs:
            a=0
            for i in tetrad.atoms:
                for j in range(3):
                    b=0
                    for k in tetrad.atoms:
                        for l in range(3):
                            H[i][j][k][l]-=\
                                tetrad.hessian[a][j][b][l]
                        b+=1
                a+=1

            # Add trimer overlap from tetrad
            for trimer in loTrimerJobs:
                match=True
                for atom in trimer.atoms:
                    if atom in tetrad.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in trimer.atoms:
                        for j in range(3):
                            b=0
                            for k in trimer.atoms:
                                for l in range(3):
                                    H[i][j][k][l]+=\
                                     trimer.hessian[a][j][b][l]
                                b+=1
                        a+=1

                    # Sub pair overlap from trimer
                    for pair in loPairJobs:
                        pairMatch=True
                        for atom in pair.atoms:
                            if atom in trimer.atoms:
                                pass
                            else:
                                pairMatch=False
                                break
                        if pairMatch:
                            a=0
                            for i in pair.atoms:
                                for j in range(3):
                                    b=0
                                    for k in pair.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]-=\
                                             pair.hessian[a][j][b][l]
                                        b+=1
                                a+=1
        
                            # Add mon overlap from pair
                            for mon in loMonJobs:
                                monMatch=True
                                for atom in mon.atoms:
                                    if atom in pair.atoms:
                                        pass
                                    else:
                                        monMatch=False
                                        break
                                if monMatch:
                                    a=0
                                    for i in mon.atoms:
                                        for j in range(3):
                                            b=0
                                            for k in mon.atoms:
                                                for l in range(3):
                                                    H[i][j][k][l]+=mon.\
                                                    hessian[a][j][b][l]
                                                b+=1
                                        a+=1

                    # Sub mon overlap from trimer
                    for mon in loMonJobs:
                        monMatch=True
                        for atom in mon.atoms:
                            if atom in trimer.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for i in mon.atoms:
                                for j in range(3):
                                    b=0
                                    for k in mon.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]-=\
                                             mon.hessian[a][j][b][l]
                                        b+=1
                                a+=1


            # Add pair overlap from tetrad
            for pair in loPairJobs:
                match=True
                for atom in pair.atoms:
                    if atom in tetrad.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in pair.atoms:
                        for j in range(3):
                            b=0
                            for k in pair.atoms:
                                for l in range(3):
                                    H[i][j][k][l]+=\
                                     pair.hessian[a][j][b][l]
                                b+=1
                        a+=1

                    # Sub mon overlap from pair
                    for mon in loMonJobs:
                        monMatch=True
                        for atom in mon.atoms:
                            if atom in pair.atoms:
                                pass
                            else:
                                monMatch=False
                                break
                        if monMatch:
                            a=0
                            for i in mon.atoms:
                                for j in range(3):
                                    b=0
                                    for k in mon.atoms:
                                        for l in range(3):
                                            H[i][j][k][l]-=\
                                             mon.hessian[a][j][b][l]
                                        b+=1
                                a+=1
            # Add mon overlap from tetrad
            for mon in loMonJobs:
                match=True
                for atom in mon.atoms:
                    if atom in tetrad.atoms:
                        pass
                    else:
                        match=False
                        break
                if match:
                    a=0
                    for i in mon.atoms:
                        for j in range(3):
                            b=0
                            for k in mon.atoms:
                                for l in range(3):
                                    H[i][j][k][l]+=\
                                        mon.hessian[a][j][b][l]
                                b+=1
                        a+=1



                
    return H
