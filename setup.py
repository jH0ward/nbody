# shared_fxns.py
import os,logging,logging.config
from subprocess import Popen,PIPE
import sys

def init(layer,input,output,msgfile):

    # Get scratch dir
    scrdir=os.path.split(input)[0]

    # Get working dir
    wrkfile=scrdir+'/workdir.txt'           # file containts dir
    fh_wrkfile=open(wrkfile,'r')
    wrkdir=fh_wrkfile.readline().strip()    # dir is on first line
    fh_wrkfile.close()

    # Save a copy of input
    os.system('cp -f ' + input + ' ' + wrkdir + '/gau2exprog.txt')

    # Clean up work directory
    extdir=wrkdir+'/ExtFiles'                # directory for external files
    indexf=scrdir+'/GauExtern.index'          # file to determine iteration #
    logf=wrkdir+'/LOGFILE'
    if os.path.exists(indexf):
        iopener=open(indexf,'r')            # find iteration number
        index=int(iopener.readline().strip())
        index+=1                            # increment iteration number
    else:
        index=1                             # or it must be the 1st iteration
        for i in (extdir,logf):
            if os.path.exists(i):
                print("\n\nERROR: Delete ",i,"or rename it")
                os.sys.exit()

    # Turn index into string for naming files
    if index>9:
        suffix=str(index)
    else:
        suffix="0"+str(index)
    # Print index to a file for the next iteration to use
    os.system('echo ' + str(index) + ' > ' + indexf)
    os.system('cp ' + indexf + ' ' + wrkdir)
    

    # Set up log file
    logging.basicConfig(filename=logf,level=logging.DEBUG)
    mylogger=logging.getLogger('init module')


    return_info=(scrdir,wrkdir,suffix)

    return return_info

def restartJob(input):
    # Get scratch dir
    scrdir=os.path.split(input)[0]

    # Get working dir
    wrkfile=scrdir+'/workdir.txt'           # file containts dir
    fh_wrkfile=open(wrkfile,'r')
    wrkdir=fh_wrkfile.readline().strip()    # dir is on first line
    fh_wrkfile.close()
    os.system('cp -f ' + input + ' ' + wrkdir + '/gau2exprog.txt')

    # Determine iteration number from ExtFiles
    ExtDir=wrkdir+'/ExtFiles/'
    indexf=scrdir+'/GauExtern.index'   # file to determine iteration #

    if os.path.exists(indexf):
        iopener=open(indexf,'r')            # find iteration number
        index=int(iopener.readline().strip())
        index+=1                            # increment iteration number
    else:
        index=1                             # or it must be the 1st iteration

    if index > 9:
        suffix=str(index)
    else:
        suffix='0'+str(index)


    os.system('echo ' + str(index) + ' > ' + indexf)
    os.system('cp ' + indexf + ' ' + wrkdir)

    # Set up log file
    logf=wrkdir+'/LOGFILE'
    logging.basicConfig(filename=logf,level=logging.DEBUG)
    mylogger=logging.getLogger('Restart')
    mylogger.info('\nJob has been restarted\n')

    return_info=(scrdir,wrkdir,suffix)

    return return_info
