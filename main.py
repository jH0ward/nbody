import sys
import os
import getpass
import datetime
from subprocess import call
import logging
import extract
import assemble
from setup import init, restartJob
from files import create
from frags import find_mons, read_mons
from COMPS import JOBS

# from COMPS_pbs import JOBS

from collections import deque

call_time = datetime.datetime.now()
# Accept input from g09
level = int(sys.argv[1])
if level == 0:
    print("Requesting external calculation of cluster")
elif level == 1:
    print("Requesting 1-body QM:QM calculation")
elif level == 2:
    print("Requesting 2-body QM:QM calculation")
elif level == 3:
    print("Requesting 3-body QM:QM calculation")
elif level == 4:
    print("Requesting 4-body QM:QM calculation")
else:
    print("Invalid option")
    sys.exit()

if sys.argv[2] == "Restart":
    print("Received restart flag")
    layer, input, output, msgfile = sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
    restart = 1
else:
    layer, input, output, msgfile = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    restart = 0

# Open input and determine the calc type
f = open(input, "r")
deriv = int(f.readline().split()[1])
f.close()

# Call initialize function
if not restart:
    return_info = init(layer, input, output, msgfile)
else:
    return_info = restartJob(input)

scrdir = return_info[0]
wrkdir = return_info[1]
suffix = return_info[2]

tfile = wrkdir + "/tflag"
f12file = wrkdir + "/f12a"
rhfile = wrkdir + "/rhf"
scalefile = wrkdir + "/scale"
corrFile = wrkdir + "/correlation"
fragFile = wrkdir + "/frags.in"
if os.path.exists(tfile):
    print("I found ", tfile)
    tflag = True
else:
    tflag = False
if os.path.exists(f12file):
    f12aFlag = True
else:
    f12aFlag = False
if os.path.exists(rhfile):
    f12aFlag = "rhf"
if os.path.exists(scalefile):
    scaleTrips = True
    f = open(scalefile, "r")
    numer = float(f.readline().split()[1])
    denom = float(f.readline().split()[1])
    trip_factor = numer / denom
    f.close()
else:
    scaleTrips = False
    trip_factor = -1.0
if os.path.exists(corrFile):
    CORR = True
else:
    CORR = False
if os.path.exists(fragFile):
    MANUAL_FRAGMENTS = True
else:
    MANUAL_FRAGMENTS = False


# Configure main.py logger
mainlogger = logging.getLogger("Iter%s" % suffix)
mainlogger.info("Iteration %d begun at %s\n\n" % (int(suffix), call_time))
mainlogger.info("Deriv = %d\n" % deriv)

# Call create function for creating input files
create(wrkdir, scrdir, suffix, level, tflag)

# Call find_mons function to map the monomers
if MANUAL_FRAGMENTS:
    mainlogger.info("Requested manual fragment specification...")
    try:
        coords, mons = read_mons(input, fragFile)
    except Exception as e:
        mainlogger.exception("Had a problem reading manual fragments: ", str(e))
        raise
else:
    mainlogger.info("Performing automatic fragment finding...")
    coords, mons = find_mons(input)

# Create monomer input
if level > 0:
    nmons = len(mons)
    hiMonJobs = []
    loMonJobs = []
    for i in range(nmons):
        hiMonJobs.append(JOBS(mons[i], coords, "hi", "mon", i))
        if tflag == False:
            loMonJobs.append(JOBS(mons[i], coords, "lo", "mon", i))

        hiMonJobs[i].makeInput(wrkdir, suffix)
        if tflag == False:
            loMonJobs[i].makeInput(wrkdir, suffix)

# Create dimer input
if level > 1:
    pairs = [i + j for i in mons for j in mons if i < j]
    npairs = len(pairs)
    hiPairJobs = []
    loPairJobs = []
    for i in range(npairs):
        hiPairJobs.append(JOBS(pairs[i], coords, "hi", "pair", i))
        if tflag == False:
            loPairJobs.append(JOBS(pairs[i], coords, "lo", "pair", i))

        hiPairJobs[-1].makeInput(wrkdir, suffix)
        if tflag == False:
            loPairJobs[-1].makeInput(wrkdir, suffix)

# Create trimer input
if level > 2:
    trimers = [i + j + k for i in mons for j in mons for k in mons if i < j < k]
    ntrimers = len(trimers)
    hiTrimerJobs = []
    loTrimerJobs = []
    for i in range(ntrimers):
        hiTrimerJobs.append(JOBS(trimers[i], coords, "hi", "trimer", i))
        if tflag == False:
            loTrimerJobs.append(JOBS(trimers[i], coords, "lo", "trimer", i))

        hiTrimerJobs[-1].makeInput(wrkdir, suffix)
        if tflag == False:
            loTrimerJobs[-1].makeInput(wrkdir, suffix)

# Create tetrad input
if level > 3:
    tetrads = [
        i + j + k + l
        for i in mons
        for j in mons
        for k in mons
        for l in mons
        if i < j < k < l
    ]
    ntetrads = len(tetrads)
    hiTetradJobs = []
    loTetradJobs = []
    for i in range(ntetrads):
        hiTetradJobs.append(JOBS(tetrads[i], coords, "hi", "tetrad", i))
        if not tflag:
            loTetradJobs.append(JOBS(tetrads[i], coords, "lo", "tetrad", i))

        hiTetradJobs[-1].makeInput(wrkdir, suffix)
        if not tflag:
            loTetradJobs[-1].makeInput(wrkdir, suffix)

# Create cluster input
natom = len(coords[0])
clusterJob = JOBS(list(range(natom)), coords, "", "cluster", "")
if not tflag:
    clusterJob.makeInput(wrkdir, suffix)
print("Colemon")


# Run everything and extract info
newdir = wrkdir + "/ExtFiles/iter_" + suffix + "/"
call(["cp", "-f", wrkdir + "/gau2exprog.txt", newdir])
# New code to determine user
user = getpass.getuser()
# process=subprocess.Popen(['qstat -u '+user+' | grep -c mem8g'],shell=True,stdout=subprocess.PIPE)
# init_jobs=int(process.communicate()[0])
# Run loMonJobs
if level > 0:
    if not tflag:
        jobQ = deque(loMonJobs)
        while len(jobQ) > 0:
            thisJob = jobQ.popleft()
            wstat = thisJob.runInput(newdir, scrdir)
    # Run hiMon Jobs
    jobQ = deque(hiMonJobs)
    while len(jobQ) > 0:
        thisJob = jobQ.popleft()
        wstat = thisJob.runInput(newdir, scrdir)

    # Run loPair Jobs
if level > 1:
    if not tflag:
        jobQ = deque(loPairJobs)
        while len(jobQ) > 0:
            thisJob = jobQ.popleft()
            wstat = thisJob.runInput(newdir, scrdir)
    #
    # Run hiPair Jobs
    wstat = 0
    jobQ = deque(hiPairJobs)
    while len(jobQ) > 0:
        thisJob = jobQ.popleft()
        wstat = thisJob.runInput(newdir, scrdir)

# Run lo trimer Jobs
if level > 2:
    wstat = 0
    if not tflag:
        jobQ = deque(loTrimerJobs)
        while len(jobQ) > 0:
            thisJob = jobQ.popleft()
            wstat = thisJob.runInput(newdir, scrdir)

    # Run hi Trimer jobs
    wstat = 0
    jobQ = deque(hiTrimerJobs)
    while len(jobQ) > 0:
        thisJob = jobQ.popleft()
        wstat = thisJob.runInput(newdir, scrdir)

# Run lo tetramer Jobs
if level > 3:
    wstat = 0
    if not tflag:
        jobQ = deque(loTetradJobs)
        while len(jobQ) > 0:
            thisJob = jobQ.popleft()
            wstat = thisJob.runInput(newdir, scrdir)

    # Run hi tetrad jobs
    jobQ = deque(hiTetradJobs)
    while len(jobQ) > 0:
        thisJob = jobQ.popleft()
        wstat = thisJob.runInput(newdir, scrdir)

# Run cluster job
if not tflag:
    clusterJob.runInput(newdir, scrdir)

if level > 0:
    if not tflag:
        for i in loMonJobs:
            #            i.runInput(newdir,scrdir)
            extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
            if deriv > 0:
                extract.grad(i, newdir)
                if deriv > 1:
                    extract.hessian(i, newdir)

if level > 1:
    if not tflag:
        for i in loPairJobs:
            #            i.runInput(newdir,scrdir)
            extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
            if deriv > 0:
                extract.grad(i, newdir)
                if deriv > 1:
                    extract.hessian(i, newdir)

if level > 2:
    if not tflag:
        for i in loTrimerJobs:
            #            i.runInput(newdir,scrdir)
            extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
            if deriv > 0:
                extract.grad(i, newdir)
                if deriv > 1:
                    extract.hessian(i, newdir)

    if level > 3:
        if not tflag:
            for i in loTetradJobs:
                print("yoo")
                print(i)
                print(i.output)
                mainlogger.info("Inside lotet\n")
                mainlogger.info("i = ", i)
                mainlogger.info("i.label = ", i.label)
                #                i.runInput(newdir,scrdir)
                extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
                if deriv > 0:
                    extract.grad(i, newdir)
                    if deriv > 1:
                        extract.hessian(i, newdir)

if level > 0:
    for i in hiMonJobs:
        extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
        if deriv > 0:
            extract.grad(i, newdir)
            if deriv > 1:
                extract.hessian(i, newdir)

if level > 1:
    for i in hiPairJobs:
        extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
        if deriv > 0:
            extract.grad(i, newdir)
            if deriv > 1:
                extract.hessian(i, newdir)

if level > 2:
    for i in hiTrimerJobs:
        extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
        if deriv > 0:
            extract.grad(i, newdir)
            if deriv > 1:
                extract.hessian(i, newdir)

    if level > 3:
        for i in hiTetradJobs:
            extract.energy(i, newdir, f12aFlag, trip_factor, CORR)
            if deriv > 0:
                extract.grad(i, newdir)
                if deriv > 1:
                    extract.hessian(i, newdir)

if not tflag:
    extract.energy(clusterJob, newdir, f12aFlag, trip_factor, CORR)
    if deriv > 0:
        extract.grad(clusterJob, newdir)
        if deriv > 1:
            extract.hessian(clusterJob, newdir)

Jobs = [clusterJob]
if level > 0:
    Jobs.extend([loMonJobs, hiMonJobs])
    if level > 1:
        Jobs.extend([loPairJobs, hiPairJobs])
        if level > 2:
            Jobs.extend([loTrimerJobs, hiTrimerJobs])
            if level > 3:
                Jobs.extend([loTetradJobs, hiTetradJobs])

# Assemble and send back
E = assemble.newenergy(Jobs, level, tflag)
if deriv > 0:
    G = assemble.newgrad(Jobs, level, tflag)
    if deriv > 1:
        H = assemble.newhess(Jobs, level, tflag)

# Send it back
f = open(output, "w")
print("%.12f," % E, "0.0, 0.0, 0.0", file=f)  ### zeroes are for dipole
if deriv > 0:
    for i in range(natom):
        print("%.14f," % G[i][0], "%.12f," % G[i][1], "%.12f" % G[i][2], file=f)

    if deriv > 1:
        for i in range(2):
            print("0.0, 0.0, 0.0", file=f)
        for i in range(3 * natom):
            print("0.0, 0.0, 0.0", file=f)
        perline = 0
        for i in range(natom):
            for j in range(3):
                for k in range(i + 1):
                    if k == i:
                        lrange = j + 1
                    else:
                        lrange = 3
                    for l in range(lrange):
                        if perline == 3:
                            print("\n", end=" ", file=f)
                            perline = 0
                        print("%.14f, " % H[i][j][k][l], end=" ", file=f)
                        perline += 1
f.close()
call(["cp", "-f", output, wrkdir + "/checkout.txt"])
