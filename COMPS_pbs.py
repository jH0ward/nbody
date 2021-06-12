# Comps.py
import numpy as np
import logging
from subprocess import check_call
import os, sys
from os import chdir


class JOBS:
    def __init__(self, atoms, coords, level, size, label):
        """
        Class of jobs to run for n-body scheme

        Methods:
            self.makeInput(wrkdir,suffix)
                -Creates input file for job
            self.runInput(wrkdir,scrdir)
                -Runs input file

        Attributes:
            self.atoms   = list of atom numbers relative to entire cluster
            self.coords  = (natom,3) array of coordinates
            self.level   = 'hi' or 'lo'
            self.size    = 'mon', 'pair',etc.
            self.label   =  level+size+<integer> (lomon0,hipair1,etc.)
            self.natom   = take a wild guess (integer)
            self.symbols = list of atomic symbols
                          ( currently not used, i think)

            NOTE: Next 2 attributes assigned during makeInput
            -------------------------------------------------
            self.exprog =  string containing name of external program
            self.cmd    =  string containing command to run exprog
        """

        self.level = level
        self.size = size
        self.label = level + size + str(label)
        self.atoms = atoms
        self.natom = len(self.atoms)
        self.symbols = []
        self.coords = np.zeros(self.natom * 3)
        self.coords.shape = (self.natom, 3)
        for i in range(self.natom):
            self.symbols.append(coords[0][atoms[i]])
            self.coords[i] = coords[1][atoms[i], :]

    def makeInput(self, wrkdir, suffix):

        JOBlogger = logging.getLogger(self.label)
        infile = wrkdir + "/ExtFiles/iter_" + suffix + "/" + self.label + ".in"

        self.input = infile
        self.base = self.input.split(".")[0]
        header = wrkdir + "/" + self.level + self.size + "_header.in"
        footer = wrkdir + "/" + self.level + self.size + "_footer.in"

        # Use command file to assign .cmd and .exprog
        cmdfile = wrkdir + "/" + self.level + self.size + "_cmd.in"
        f_cmd = open(cmdfile, "r")
        self.exprog = f_cmd.readline().strip()
        self.cmd = f_cmd.readline().strip()
        f_cmd.close()

        # If this is a restart, get the hell out
        if os.path.exists(infile):
            JOBlogger.info("\nSkipping %s\n" % infile)
            return 0

        # Input for writing and head & foot for reading
        f_in = open(infile, "w")
        f_head = open(header, "r")
        f_foot = open(footer, "r")

        # Drop header
        for line in f_head:
            print(line.rstrip(), file=f_in)
        f_head.close()

        # Drop coordinates
        for i in range(self.natom):
            # Print atomic symbol 1st
            print(self.symbols[i], end=" ", file=f_in)

            # NOTE: mpqc needs extra [
            if self.exprog == "mpqc":
                print("[", end=" ", file=f_in)

            # Print x,y,z for atom i
            for j in self.coords[i]:
                print("%.12f" % j, end=" ", file=f_in)

            # NOTE: mpqc needs extra ]
            if self.exprog == "mpqc":
                print("]", end=" ", file=f_in)

            # Move to next line
            print("\n", end=" ", file=f_in)

        # Drop footer
        for line in f_foot:
            print(line, end=" ", file=f_in)
        f_foot.close()
        f_in.close()

        JOBlogger.info("\nNew Input created --> %s\n" % self.input)

        return 0

    def runInput(self, wrkdir, scrdir):
        JOBlogger = logging.getLogger(self.label)
        base = self.base
        self.output = base + ".out"
        newdir = self.base.split("/")[:-1]

        # if .out exists don't fuck with it
        if os.path.exists(self.output):
            ## Set g09 output to formatted chkpt file
            # if self.exprog=='g09':
            #    self.output=base+'.fchk'
            JOBlogger.info("\nSkipping %s", self.output)
            return 0

        chdir(wrkdir)

        # /qc/bin/bin_cfour_v1.sh syntax
        if self.exprog == "cfour":
            # check_call(['cp','-f',self.input,'ZMAT'])
            c4file = self.input.split("/")[-1]
            cmd = self.cmd.split()[0]
            nproc = self.cmd.split()[1]
            # check_call([cmd,nproc])
            syntax = cmd + " " + c4file + " " + nproc
            # check_call([cmd,c4file,nproc])

        # pqs4.sh syntax
        elif self.exprog == "pqs":
            cmd = self.cmd.split()[0]
            nproc = self.cmd.split()[1]
            # check_call([cmd,nproc,self.input])
            pqs_f = self.input.split("/")[-1]
            syntax = cmd + " " + nproc + " " + pqs_f

        # mpqc.sh and g09 syntax
        else:
            g09file = self.input.split("/")[-1]
            cmd = self.cmd.split()[0]
            # check_call([self.cmd,g09file])
            # check_call([cmd,g09file])
            syntax = cmd + " " + g09file

        # Define pbs_file and open for writing
        pbs_fname = base + ".pbs"
        pbs_f = open(pbs_fname, "w")

        # TODO Stop printing full input file name to pbs file, this is wrong
        print(
            """\
#!/bin/bash
#
#PBS -N """
            + self.label
            + """
#PBS -l mem=7200mb
#PBS -l ncpus=6
#PBS -j oe
#PBS -W umask=022
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR

"""
            + syntax
            + """

""",
            file=pbs_f,
        )
        pbs_f.close()

        # Submit that job
        check_call(["qsub", pbs_fname])
        # clean up cfour output
        if self.exprog == "cfour":
            # check_call(['mv','./JOBDUMP',base+'.JOBDUMP'])
            # check_call(['mv','./output.dat',base+'.out'])
            check_call(["rm", "-f", "./MOLDEN"])
            check_call(["rm", "-f", "./DIPOL"])
            # if os.path.exists('./FCMFINAL'):
            #    check_call(['mv','./FCMFINAL',base+'.FCMFINAL'])
            if os.path.exists("./FCMINT"):
                check_call(["rm", "-f", "./FCMINT"])

        # clean up pqs output
        # if self.exprog=='pqs':
        # TODO move this section to code that runs after data collection
        # because this shit doesn't even exist yet
        # print "CaN YOU THINK OF AN ANIMAL THAT TALKS..BESIDES A PERSON"
        # print "Coleman trying to delete", base+'.in.out'
        # print "Coleman trying to delete", base+'.basis'
        # print "Coleman trying to delete", base+'.mos'
        # check_call(['rm','-f',base+'.in.out'])
        # check_call(['rm','-f',base+'.basis'])
        # check_call(['rm','-f',base+'.basis2'])
        ##check_call(['rm','-f',base+'.control'])
        # check_call(['rm','-f',base+'.in.log'])
        # check_call(['rm','-f',base+'.in.messages'])
        # check_call(['rm','-f',base+'.mos'])
        # check_call(['rm','-f',base+'.log'])

        # clean up mpqc output
        if os.path.exists(base + ".wfn"):
            check_call(["gzip", "-f", base + ".wfn"])

        # for g09, create .chk
        # if self.exprog=='g09':
        #     tscript='/home/jumper/tschumpr/bin/g09fchk.sh'
        #     short=base.split('/')[-1]
        #     chk_fz=short+'.chk.gz'
        #     chk_f=short+'.chk'
        #     #check_call(['gunzip','ext.chk.gz'])         # unzip
        #     check_call(['gunzip',chk_fz])         # unzip
        #     #check_call([tscript,'ext.chk'])             # make .fchk
        #     check_call([tscript,chk_f])# make .fchk
        # check_call(['rm','-f','ext.chk'])           # rm .chk
        # check_call(['mv','ext.fchk',base+'.fchk'])   # rename .fchk

        chdir(scrdir)
        JOBlogger.info("\nNew Output created --> %s\n" % self.output)
        return 1
