# extract.py
from collections import deque
import os
from os import chdir
import numpy as np
import rotate
import hessian_rotate

np.set_printoptions(precision=12)
import logging
from subprocess import check_call


def energy(job, newdir, f12aFlag, scaleTrips, just_corr):
    chdir(newdir)

    exlogger = logging.getLogger(job.label)
    exlogger.info("Extacting from this job: %s" % job.label)

    output = job.output
    natom = job.natom

    # MPQC energy
    if job.exprog == "mpqc":
        f = open(output, "r")
        while True:
            line = f.readline().strip()
            if f12aFlag == "rhf":
                if "RHF energy [au]:" in line:
                    job.energy = float(line.split(":")[1])
                    exlogger.info("\nEnergy = %3.11f\n" % job.energy)
                    break
            else:
                if "Value of the MolecularEnergy:" in line:
                    job.energy = float(line.split(":")[1])
                    exlogger.info("\nEnergy = %3.11f\n" % job.energy)
                    break

    # CFOUR energy
    elif job.exprog == "cfour":
        job.output = output.split(".")[0] + ".JOBDUMP"
        output = job.output
        f = open(output, "r")
        job.energy = float(f.readline().split()[1])
        exlogger.info("\nEnergy = %3.11f\n" % job.energy)

    # PQS energy
    elif job.exprog == "pqs":
        # need to change output file from grabbing energy to .control
        output = job.base + ".control"
        f = open(output, "r")
        while True:
            line = f.readline().strip()
            # if 'Total Energy =' in line:
            if "$energy" in line:
                # job.energy=float(line.split()[3])
                job.energy = float(line.split()[1])
                exlogger.info("\nEnergy = %3.11f\n" % job.energy)
                break
        base = job.base
        # Now it's time to clean up pqs unncecessary bullshit
        check_call(["rm", "-f", base + ".in.out"])
        check_call(["rm", "-f", base + ".basis"])
        check_call(["rm", "-f", base + ".basis2"])
        check_call(["rm", "-f", base + ".in.log"])
        check_call(["rm", "-f", base + ".in.messages"])
        check_call(["rm", "-f", base + ".mos"])
        check_call(["rm", "-f", base + ".log"])
        # check_call(['rm','-f',base+'.control'])
        check_call(["rm", "-f", base + ".sym"])

    elif job.exprog == "molpro":
        output = job.base + ".out"
        print("output = ", output)
        f = open(output, "r")
        if scaleTrips > -1.0:
            if f12aFlag == True:
                searchTimes = 1
            else:
                searchTimes = 2
            for i in range(searchTimes):
                while True:
                    line = f.readline()
                    if line == "":
                        break
                    if "Triples (T) contribution" in line:
                        trip = float(line.split()[-1])
                        trip_scaled = scaleTrips * trip
                        break
                f.readline()
                f.readline()
                restofE = float(f.readline().split()[-1])
            job.energy = trip_scaled + restofE
            exlogger.info("\nEnergy = %3.11f\n" % job.energy)
            f.close()
            return 0

        if just_corr:
            if f12aFlag == True:
                searchTimes = 1
            else:
                searchTimes = 2
            for i in range(searchTimes):
                while True:
                    line = f.readline()
                    if line == "":
                        break
                    if "Total correlation energy" in line:
                        corr = float(line.split()[-1])
                        break
            job.energy = corr
            exlogger.info("\nCorr. energy = %3.11f\n" % job.energy)
            f.close()
            return 0

        while True:
            line = f.readline()
            if line == "":
                break
            line = line.rstrip()
            if "!" in line:
                eline = line
                if f12aFlag == True:
                    if "!CCSD(T)-F12a total energy" in eline:
                        break
                if f12aFlag == "rhf":
                    if "!RHF STATE 1.1 Energy" in eline:
                        break
        exlogger.info("Energy grabbed from  line =\n%s\n" % eline)
        ary = eline.split()
        job.energy = float(ary[-1])
        exlogger.info("\nEnergy = %3.11f\n" % job.energy)

    # G09 energy
    elif job.exprog == "g09":
        tscript = "/home/jumper/tschumpr/bin/g09fchk.sh"
        base = job.base
        if not os.path.exists(base + ".fchk"):
            short = base.split("/")[-1]
            chk_fz = short + ".chk.gz"
            chk_f = short + ".chk"
            check_call(["gunzip", chk_fz])  # unzip
            check_call([tscript, chk_f])  # make .fchk
            check_call(["rm", "-f", chk_f])  # rm .chk
        chkf_name = base + ".fchk"
        f = open(chkf_name, "r")
        while True:
            line = f.readline()
            if "Total Energy" in line:
                arry = line.split()
                job.energy = float(arry[3])
                exlogger.info("\nEnergy = %3.14f\n" % job.energy)
                break
    f.close()
    return 0


def grad(job, newdir):
    exlogger = logging.getLogger(job.label)

    output = job.output
    natom = job.natom

    grad = np.zeros(job.natom * 3)

    # MPQC Gradient
    if job.exprog == "mpqc":
        f = open(output, "r")
        while True:
            line = f.readline().strip()
            if "Gradient of the MolecularEnergy:" in line:
                break

        for i in range(natom * 3):
            line = f.readline().strip().split()
            grad[i] = float(line[1])

        grad.shape = (natom, 3)
        job.grad = grad
        exlogger.info("\nGradient from %s\n%s" % (job.exprog, job.grad))

    # CFOUR Gradient
    elif job.exprog == "cfour":
        grad.shape = (natom, 3)
        output = job.output
        f = open(output, "r")
        # Skip to the atom map
        skip = 4 + 2 * job.natom
        for i in range(skip):
            f.readline()
        arry = f.readline().split()
        atmap = []
        for i in range(job.natom):
            atmap.append(int(arry[i]) - 1)
        f.close()
        # Now go back in and get the stuff
        f = open(output, "r")
        skip = 1 + job.natom
        for i in range(skip):
            f.readline()
        for i in range(job.natom):
            arry = f.readline().split()
            for j in range(3):
                grad[atmap[i]][j] = float(arry[j])
        grad.shape = (natom, 3)
        job.grad = grad
        exlogger.info("\nGradient from %s\n%s" % (job.exprog, job.grad))
        f.close()

    # PQS Gradient
    elif job.exprog == "pqs":
        grad.shape = (natom, 3)
        output = job.base + ".grad"
        f = open(output, "r")
        f.readline()

        for i in range(natom):
            arry = f.readline().split()
            for j in range(3):
                grad[i][j] = -1 * float(arry[j])

        job.grad = grad
        exlogger.info("\nGradient from %s\n%s" % (job.exprog, job.grad))
        f.close()

    # G09 gradient
    elif job.exprog == "g09":
        output = job.output
        grad.shape = (natom, 3)
        f = open(output, "r")
        exlogger.info("\nOutput is %s\n" % job.output)
        while True:
            line = f.readline()
            if "First derivatives punched:" in line:
                break

        for i in range(natom):
            arry = f.readline().split()
            for j in range(3):
                arry[j] = arry[j].replace("D", "E")
                grad[i][j] = float(arry[j])

        job.grad = grad
        exlogger.info("\nGradient from %s\n%s" % (job.exprog, job.grad))
        f.close()

    if not job.exprog == "g09":
        coords(job, newdir)
        ## Do rotation function
        job.rgrad = rotate.matrix(job)
        exlogger.info("\nRotated gradient:\n%s" % job.rgrad)

    else:
        job.rgrad = job.grad
        # Don't rotate g09 derivatives

    return 0


def hessian(job, newdir):
    exlogger = logging.getLogger(job.label)
    output = job.output
    natom = job.natom
    hessian = np.zeros(natom ** 2 * 9)

    # CFOUR Hessian
    if job.exprog == "cfour":
        output = job.base + ".FCMFINAL"
        f = open(output, "r")
        f.readline()
        hessian.shape = (natom, 3, natom, 3)
        for i in range(natom):
            for j in range(3):
                for k in range(natom):
                    arry = f.readline().split()
                    for l in range(3):
                        hessian[i][j][k][l] = arry[l]

        job.hessian = hessian
        tempHess = np.zeros(natom * natom * 3 * 3)
        tempHess.shape = (natom, natom, 3, 3)
        job.rhess = np.zeros(natom * natom * 3 * 3)
        job.rhess.shape = (natom, natom, 3, 3)

        # Fill temporary hessian with more convenient shape
        for i in range(natom):
            for j in range(natom):
                for k in range(3):
                    for l in range(3):
                        tempHess[i][j][k][l] = hessian[i][k][j][l]
                job.rhess[i][j] = hessian_rotate.matrix(job, tempHess[i][j])

        # Translate temporary hessian back into natom x 3 x natom x 3
        for i in range(natom):
            for j in range(3):
                for k in range(natom):
                    for l in range(3):
                        job.hessian[i][j][k][l] = job.rhess[i][k][j][l]

        f.close()
        exlogger.info("\nHessian from %s\n%s" % (job.exprog, job.hessian))

        return 0

    # g09 Hessian
    if job.exprog == "g09":
        f = open(output, "r")
        hessian.shape = (natom, 3, natom, 3)
        while True:
            line = f.readline()
            # if 'Cartesian Force Constants' in line:
            #    break
            if "Second derivatives punched:" in line:
                break

        line = f.readline()
        line = line.replace("D", "E")
        # arry=f.readline().split()
        arry = line.split()
        queue = deque(arry)
        for i in range(natom):
            for j in range(3):
                for k in range(i + 1):
                    if k == i:
                        lrange = j + 1
                    else:
                        lrange = 3
                    for l in range(lrange):
                        if len(queue) == 0:
                            line = f.readline()
                            # arry=f.readline().split()
                            if line == "\n":
                                break
                            line = line.replace("D", "E")
                            arry = line.split()
                            queue = deque(arry)
                        # print "COLEMAN QUEUE!!"
                        # print queue.popleft()
                        hessian[i][j][k][l] = queue.popleft()
                        hessian[k][l][i][j] = hessian[i][j][k][l]
                        # if arry[0]=='Dipole':
                        #    break

        job.hessian = hessian

        # Don't rotate g09 derivatives

        # tempHess=np.zeros(natom*natom*3*3)
        # tempHess.shape=(natom,natom,3,3)
        # job.rhess=np.zeros(natom*natom*3*3)
        # job.rhess.shape=(natom,natom,3,3)

        # Fill temporary hessian with more convenient shape
        # for i in range(natom):
        #    for j in range(natom):
        #        for k in range(3):
        #            for l in range(3):
        #                tempHess[i][j][k][l]=hessian[i][k][j][l]
        #        job.rhess[i][j]=hessian_rotate.matrix(job,tempHess[i][j])

        # Translate temporary hessian back into natom x 3 x natom x 3
        # for i in range(natom):
        #    for j in range(3):
        #        for k in range(natom):
        #            for l in range(3):
        #                job.hessian[i][j][k][l]=job.rhess[i][k][j][l]

        f.close()
        exlogger.info("\nHessian from %s\n%s" % (job.exprog, job.hessian))

        return 0


def coords(job, newdir):
    exlogger = logging.getLogger(job.label)
    output = job.output
    natom = job.natom

    exyz = np.zeros(3 * natom)

    # MPQC coordinates
    if job.exprog == "mpqc":
        exyz.shape = (natom, 3)
        f = open(output, "r")
        while True:
            line = f.readline().strip()
            if "geometry" in line:
                break
        for i in range(natom):
            arry = f.readline().split()
            for j in range(3):
                arry[j + 3] = arry[j + 3].strip("]")
                exyz[i][j] = arry[j + 3]

        exlogger.info("\nCoordinates from %s\n%s\n" % (job.exprog, exyz))

    # CFOUR coordinates
    elif job.exprog == "cfour":
        exyz.shape = (natom, 3)
        f = open(output, "r")
        f.readline()

        for i in range(natom):
            arry = f.readline().split()
            for j in range(3):
                exyz[i][j] = float(arry[j])

        exlogger.info("\nCoordinates from %s\n%s\n" % (job.exprog, exyz))

    # PQS coordinates
    elif job.exprog == "pqs":
        exyz.shape = (natom, 3)
        output = job.base + ".coord"
        f = open(output, "r")
        f.readline()

        for i in range(natom):
            arry = f.readline().split()
            for j in range(3):
                exyz[i][j] = float(arry[j + 1])

        exlogger.info("\nCoordinates from %s\n%s\n" % (job.exprog, exyz))

    elif job.exprog == "g09":
        output = job.output
        f = open(output, "r")
        while True:
            line = f.readline()
            if "Current cartesian coordinates" in line:
                break

        j = 0  # placeholder for coords container
        while True:
            arry = f.readline().split()
            if arry[0] == "Force":
                break
            for i in range(len(arry)):
                exyz[j] = float(arry[i])
                j += 1

        exyz.shape = (natom, 3)

    job.exyz = exyz  # job.exyz will be external xyz coordinates (after rotation)
    f.close()
    return 0
