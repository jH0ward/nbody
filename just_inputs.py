#!/usr/bin/python
import sys,os
import datetime
import time
from subprocess import call
import logging,extract,assemble,subprocess
from setup import init,restartJob
from files import create
from frags import find_mons
from COMPS import JOBS
import glob

from collections import deque

call_time=datetime.datetime.now()
# Accept input from g09
print "Here comes the args"
print sys.argv
level=int(sys.argv[1])
if level==1:
    print "Requesting external calculation on entire cluster"
elif level==2:
    print "Requesting 2-body QM:QM calculation"
elif level==3:
    print "Requesting 3-body QM:QM calculation"
else:
    print "Invalid option"
    sys.exit()

if sys.argv[2]=='Restart':
    print "Received restart flag"
    layer,input,output,msgfile=sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]
    restart=1
else:
    layer,input,output,msgfile=sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
    restart=0

# Open input and determine the calc type
call(["cp", "-f", input, '/home/jumper/jchoward/python/charge/use_g09'])
call(["cp", "-f", sys.argv[6], '/home/jumper/jchoward/python/charge/use_g09'])
call(["cp", "-f", sys.argv[7], '/home/jumper/jchoward/python/charge/use_g09'])
call(["cp", "-f", sys.argv[8], '/home/jumper/jchoward/python/charge/use_g09'])
call(["cp", "-f", sys.argv[7], '/home/jumper/jchoward/python/charge/use_g09'])
#call(["cp", "-f", '\*inp', '/home/jumper/jchoward/python/charge/use_g09'])
#call(["cp -f *inp /home/jumper/jchoward/python/charge/use_g09"])
proc=subprocess.Popen('cp *.inp /home/jumper/jchoward/python/charge/use_g09',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
out,err=proc.communicate()
print "Here comes out and error"
print out
print err

print "Contents of scratch dir:"
print glob.glob('./*')


f=open(input,'r')
deriv=int(f.readline().split()[1])
f.close()

# Call initialize function
if not restart:
    return_info = init(layer,input,output,msgfile)
else:
    return_info = restartJob(input)

scrdir=return_info[0]
wrkdir=return_info[1]
suffix=return_info[2]

# Configure main.py logger
mainlogger=logging.getLogger('Iter%s' %suffix)
mainlogger.info('Iteration %d begun at %s\n\n'%(int(suffix),call_time))
mainlogger.info("Deriv = %d\n" %deriv)

tfile=wrkdir+'/tflag'
if os.path.exists(tfile):
    tflag=True
else:
    tflag=False
# Call create function for creating input files
create(wrkdir,scrdir,suffix,level,tflag)

# Call find_mons function to map the monomers
coords,mons=find_mons(input)

# Create monomer input
if level>1:
    nmons=len(mons)
    hiMonJobs=[]
    loMonJobs=[]
    for i in range(nmons):
        hiMonJobs.append(JOBS(mons[i],coords,'hi','mon',i))
        loMonJobs.append(JOBS(mons[i],coords,'lo','mon',i))

        hiMonJobs[i].makeInput(wrkdir,suffix)
        loMonJobs[i].makeInput(wrkdir,suffix)

# Create dimer input
    pairs=[ i + j for i in mons for j in mons if i<j ]
    npairs=len(pairs)
    hiPairJobs=[]
    loPairJobs=[]
    for i in range(npairs):
        hiPairJobs.append(JOBS(pairs[i],coords,'hi','pair',i))
        loPairJobs.append(JOBS(pairs[i],coords,'lo','pair',i))

        hiPairJobs[-1].makeInput(wrkdir,suffix)
        loPairJobs[-1].makeInput(wrkdir,suffix)

# Create trimer input
if level>2:
    trimers=[ i + j + k for i in mons for j in mons for k in mons if i<j<k ]
    ntrimers=len(trimers)
    hiTrimerJobs=[]
    loTrimerJobs=[]
    for i in range(ntrimers):
        hiTrimerJobs.append(JOBS(trimers[i],coords,'hi','trimer',i))
        loTrimerJobs.append(JOBS(trimers[i],coords,'lo','trimer',i))

        hiTrimerJobs[-1].makeInput(wrkdir,suffix)
        loTrimerJobs[-1].makeInput(wrkdir,suffix)

# Create cluster input
natom=len(coords[0])
clusterJob=JOBS(range(natom),coords,'','cluster','')
clusterJob.makeInput(wrkdir,suffix)

# Coleman only edit from main.py
sys.exit()

# Run everything and extract info
newdir=wrkdir+'/ExtFiles/iter_'+suffix+'/'
process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
    shell=True,stdout=subprocess.PIPE)
init_jobs=int(process.communicate()[0])
# Run loMonJobs
#if level>1:
#    jobQ=deque(loMonJobs)
#    while len(jobQ)>0:
#        #process=subprocess.Popen(['qstat -u jchoward | grep -c mem4g'],
#        #    shell=True,stdout=subprocess.PIPE)
#        process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#            shell=True,stdout=subprocess.PIPE)
#        nJobs=int(process.communicate()[0])
#        #if nJobs>5:
#        if nJobs>init_jobs+5:
#            time.sleep(5)
##        else:
#            thisJob=jobQ.popleft()
#            wstat=thisJob.runInput(newdir,scrdir)
#            if wstat:
#                time.sleep(2)
#            else:
#                mainlogger.info("Not waiting on lo mon jobs\n")
#
#    wstat=0
#    # Run hiMon Jobs
#    jobQ=deque(hiMonJobs)
#    while len(jobQ)>0:
#        #process=subprocess.Popen(['qstat -u jchoward | grep -c mem4g'],
#        #    shell=True,stdout=subprocess.PIPE)
#        process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#            shell=True,stdout=subprocess.PIPE)
#        nJobs=int(process.communicate()[0])
#        if nJobs>init_jobs+5:
#            time.sleep(5)
#        else:
#            thisJob=jobQ.popleft()
#            wstat=thisJob.runInput(newdir,scrdir)
#            if wstat:
#                time.sleep(2)
#            else:
#                mainlogger.info("Not waiting on hi mon jobs\n")

    # Run loPair Jobs
#    wstat=0
#    jobQ=deque(loPairJobs)
#    while len(jobQ)>0:
#        process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#            shell=True,stdout=subprocess.PIPE)
#        nJobs=int(process.communicate()[0])
#        if nJobs>init_jobs+5:
#            time.sleep(5)
#        else:
#            thisJob=jobQ.popleft()
#            wstat=thisJob.runInput(newdir,scrdir)
#            if wstat:
#                time.sleep(2)
#            else:
#                mainlogger.info("Not waiting on lo pair jobs\n")
#
    # Run hiPair Jobs
#    wstat=0
#    jobQ=deque(hiPairJobs)
#    while len(jobQ)>0:
#        process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#            shell=True,stdout=subprocess.PIPE)
#        nJobs=int(process.communicate()[0])
#        if nJobs>init_jobs+5:
#            time.sleep(5)
#        else:
#            thisJob=jobQ.popleft()
##            wstat=thisJob.runInput(newdir,scrdir)
#            if wstat:
#                time.sleep(2)
#            else:
#                mainlogger.info("Not waiting on hi pair jobs\n")

# Run lo trimer Jobs
#if level>2:
#    wstat=0
#    jobQ=deque(loTrimerJobs)
#    while len(jobQ)>0:
#        process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#            shell=True,stdout=subprocess.PIPE)
#        nJobs=int(process.communicate()[0])
#        if nJobs>init_jobs+5:
#            time.sleep(5)
#        else:
#            thisJob=jobQ.popleft()
#            wstat=thisJob.runInput(newdir,scrdir)
#            if wstat:
#                time.sleep(2)
#            else:
#                mainlogger.info("Not waiting on lo trimer jobs\n")

    # Run hi Trimer jobs
#    wstat=0
#    jobQ=deque(hiTrimerJobs)
#    while len(jobQ)>0:
#        process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#            shell=True,stdout=subprocess.PIPE)
#        nJobs=int(process.communicate()[0])
#        if nJobs>init_jobs+5:
#            time.sleep(5)
#        else:
#            thisJob=jobQ.popleft()
#            wstat=thisJob.runInput(newdir,scrdir)
#            if wstat:
#                time.sleep(2)
#            else:
#                mainlogger.info("Not waiting on hi trimer jobs\n")
#
#    wstat=0
#
## Run cluster job
#clusterJob.runInput(newdir,scrdir)

#wait=1
#while wait:
#    process=subprocess.Popen(['qstat -u jchoward | grep -c mem8g'],
#        shell=True,stdout=subprocess.PIPE)
#    nJobs=int(process.communicate()[0])
#    if nJobs>init_jobs:
#        time.sleep(10)
#    else:
#        wait=0
#

if level>1:
    for i in loMonJobs:
        i.runInput(newdir,scrdir)
        extract.energy(i,newdir)
        if deriv>0:
            extract.grad(i,newdir)
            if deriv>1:
                extract.hessian(i,newdir)

    for i in loPairJobs:
        i.runInput(newdir,scrdir)
        extract.energy(i,newdir)
        if deriv>0:
            extract.grad(i,newdir)
            if deriv>1:
                extract.hessian(i,newdir)

if level>2:
    for i in loTrimerJobs:
        i.runInput(newdir,scrdir)
        extract.energy(i,newdir)
        if deriv>0:
            extract.grad(i,newdir)
            if deriv>1:
                extract.hessian(i,newdir)

if level>1:
    for i in hiMonJobs:
        i.runInput(newdir,scrdir)
        extract.energy(i,newdir)
        if deriv>0:
            extract.grad(i,newdir)
            if deriv>1:
                extract.hessian(i,newdir)

    for i in hiPairJobs:
        i.runInput(newdir,scrdir)
        extract.energy(i,newdir)
        if deriv>0:
            extract.grad(i,newdir)
            if deriv>1:
                extract.hessian(i,newdir)

if level>2:
    for i in hiTrimerJobs:
        i.runInput(newdir,scrdir)
        extract.energy(i,newdir)
        if deriv>0:
            extract.grad(i,newdir)
            if deriv>1:
                extract.hessian(i,newdir)

clusterJob.runInput(newdir,scrdir)
extract.energy(clusterJob,newdir)
if deriv>0:
    extract.grad(clusterJob,newdir)
    if deriv>1:
        extract.hessian(clusterJob,newdir)

Jobs=[clusterJob]
if level>1:
    Jobs.extend([loMonJobs,hiMonJobs,loPairJobs,hiPairJobs])
    if level>2:
        Jobs.extend([loTrimerJobs,hiTrimerJobs])

# Assemble and send back
E=assemble.energy(Jobs,level)
if deriv>0:
    G=assemble.gradient(Jobs,level)
    if deriv>1:
        H=assemble.hessian(Jobs,level)

 ## Send it back
f= open(output, 'w')
print >> f, "%.12f," % E, "0.0, 0.0, 0.0"  ### zeroes are for dipole
if deriv>0:
    for i in range(natom):
        print >> f, "%.14f," % G[i][0], "%.12f," % G[i][1],"%.12f" % G[i][2]

    if deriv>1:
        for i in range(2):
            print >> f, "0.0, 0.0, 0.0"
        for i in range(3*natom):
            print >> f, "0.0, 0.0, 0.0"
        perline=0
        for i in range(natom):
            for j in range(3):
                for k in range(i+1):
                    if k==i:
                        lrange=j+1
                    else:
                        lrange=3
                    for l in range(lrange):
                        if perline==3:
                            print >> f, "\n",
                            perline=0
                        print >> f, "%.14f, "%H[i][j][k][l],
                        perline+=1
call(["cp", "-f", output, "/home/jumper/jchoward/checkout.txt"])
f.close()
