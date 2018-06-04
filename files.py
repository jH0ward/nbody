#files.py
import os,sys,logging
from subprocess import call

def create(wrkdir,scrdir,suffix,level,tflag):

    createlogger=logging.getLogger('Files')
    if level==1:
        levels = ['lo','hi']
        if tflag==True:
            levels=['hi']
        types = ['mon']
        templates = ['_header.in', '_footer.in']

    if level>1:
        levels = ['lo','hi']
        if tflag==True:
            levels=['hi']
        types = ['mon','pair']
        templates = ['_header.in', '_footer.in']
    else:
        levels=[]
        types=[]
        templates=[]
    if level>2:
        types.append('trimer')
    if level>3:
        types.append('tetrad')



    # Make sure headers and footers exist
    head_and_foot = [wrkdir+'/'+ i + j + k for i in levels \
                                           for j in types  \
                                           for k in templates ]
    # Check for monomer and pair header and footer
    missing = [i for i in head_and_foot if not os.path.exists(i) ]
    if missing:
        print "You are missing ", missing
        print "Exiting now..."
        createlogger.critical('\nMissing file %s\nExiting....\n' %missing)
        sys.exit()

    # Check for cluster header and footer
    if tflag==False:
        cluster_files = [wrkdir+'/cluster'+i for i in templates]

        missing = [i for i in cluster_files if not os.path.exists(i) ]
        if missing:
            print "You are missing ",missing
            print "Exiting now..."
            createlogger.critical('\nMissing file %s\nExiting....' %missing)
            sys.exit()

    # Check for cmd files
    cmd_files = [wrkdir+'/'+ i + j +'_cmd.in' for i in levels \
                                              for j in types ]
    if tflag==False:
        cmd_files.append(wrkdir+'/cluster_cmd.in')
    
    missing = [i for i in cmd_files if not os.path.exists(i)]
    if missing:
        print "You are missing ", missing
        createlogger.critical('\nMissing file %s\nExiting....\n' %missing)
        sys.exit()
        print "Exiting now..."
    newdir = wrkdir+'/ExtFiles/iter_'+suffix
    call(['mkdir', '-p', newdir ])
    call(['cp','-f',wrkdir+'/GENBAS',newdir])

    if os.path.exists(wrkdir+'/ECPDATA'):
        call(['cp','-f',wrkdir+'/ECPDATA',newdir])

    # copy hostfile
    if os.path.exists(wrkdir+'/hostfile'):
        call(['cp','-f',wrkdir+'/hostfile',newdir])
        createlogger.info('\nCopied hostfile' )

    return 0
