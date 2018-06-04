# rotate.py
import numpy
from numpy import linalg,dot,mat,zeros,cross
import logging,math,sys
from math import fabs

def matrix(job,mu=False,dip=False):
    rotlogger=logging.getLogger(job.label)
    output=job.output
    natom=job.natom
    xyz_gau=job.coords
    xyz_ext=job.exyz
    grad_ext=job.grad

    small = 1e-8
    numpy.set_printoptions(precision=15)

    Txyz_gau = numpy.zeros(natom*3)
    Txyz_gau.shape = (natom,3)


    Txyz_ext = numpy.zeros(natom*3)
    Txyz_ext.shape = (natom,3)


### Check for translation / rotation
    max=numpy.zeros(3)
    diff=numpy.zeros(3)
    last=numpy.zeros(3)
    do_i_rot = 'no'
    for i in range (1, natom):
        for j in range (3):
            diff[j] = fabs((xyz_ext[i][j] - \
               xyz_gau[i][j])-(xyz_ext[i-1][j]-xyz_gau[i-1][j]))
            if diff[j] > small:
                do_i_rot = "yeah"
             
    
    
    if do_i_rot == "yeah":
        rotlogger.warning("\nRotation Detected\n")
        ### Translate both geometries to put atom 1 at the origin
        for i in range (natom):
            for j in range (3):
                Txyz_ext[i][j] = xyz_ext[i][j]-xyz_ext[0][j]
                Txyz_gau[i][j] = xyz_gau[i][j]-xyz_gau[0][j]
           
        ### Make unit vector from atom1 to atom 2
        gau_v = []
        ext_v = []
        for i in range(3):
            gau_v.append(Txyz_gau[1][i]/linalg.norm(Txyz_gau[1]))
            ext_v.append(Txyz_ext[1][i]/linalg.norm(Txyz_ext[1]))
        
        ### Dot the two vectors
        costh = dot(gau_v, ext_v)
        ### Cross and normalize
        w=[]
        w = numpy.cross(gau_v, ext_v)
        w_u=[]
        w_r=[]
        ### Unit vector of rotation axis
        for i in range(3):
            if linalg.norm(w)!=0:
                w_u.append(w[i]/linalg.norm(w))
            else:
                w_u.append(w[i])
            w_r.append(-1*w[i]/linalg.norm(w)) #antiparallel
        sinth = linalg.norm(w)
        #### 1st part of rotation matrix
        tmp1rot = RodRot(costh,sinth,w_u)
        ### Test back rotation for atom 2
        ext2 = numpy.zeros(3)
        ext2 = Txyz_ext[1]
        ### Default rotation
        M_ext2 = mat(ext2.copy())
        M_tmp1rot = mat(tmp1rot.copy())
        tmp1 = M_ext2 * M_tmp1rot
        ### Check difference
        mdiff = 0.0
        #ab = math.fabs(neg1)
        diff=zeros(3)
        a_tmp1 = numpy.asarray(tmp1)
        for i in range(3):
            diff[i] = (a_tmp1[0][i] - Txyz_gau[1][i])
            if fabs(diff[i]) > mdiff:
                mdiff = fabs(diff[i])
        
        rotlogger.info("Max diff on default rotation = %3.11f\n" %mdiff)
        #if mdiff < small:
        #    rotlogger.info("\nExternal coordinates fixed by default rotation\n")
           #print "Default rotation works"
        
        if not mdiff < small:
            tmp1rot = RodRot(costh,sinth,w_r) #Reverse rot
            M_tmp1rot = mat(tmp1rot.copy())
            tmp2 = M_ext2 * M_tmp1rot
            mdiff = 0.0
            diff=zeros(3)
            a_tmp2 = numpy.asarray(tmp2)
            for i in range(3):
               diff[i] = (a_tmp2[0][i] - Txyz_gau[1][i])
               if fabs(diff[i]) > mdiff:
                  mdiff = fabs(diff[i])
           
           #print "Max diff on reverse rotation is ", mdiff
            #if mdiff > small:
            #    rotlogger.warning("\nRotation not fixed yet\n")

        
        ## Perform 1st rotation
        M_Txyz_ext = mat(Txyz_ext.copy())
        M_tmp1rot = mat(tmp1rot.copy())
        tmp_xyz = M_Txyz_ext * M_tmp1rot
        a_tmp_xyz = numpy.asarray(tmp_xyz)
        for i in range (2):
            for j in range (3):
                diff[j] = a_tmp_xyz[i][j] - Txyz_gau[i][j]
        
        ### Find 2nd part of rotation matrix
        ### using coords of atom 2 as axis
        w=zeros(3)
        w[0] = a_tmp_xyz[1][0]
        w[1] = a_tmp_xyz[1][1]
        w[2] = a_tmp_xyz[1][2]
        w_u=[]
        w_r=[]
        
        for i in range(3):
            w_u.append(w[i]/linalg.norm(w))
            w_r.append(-1*w[i]/linalg.norm(w)) #antiparallel
        sinth = linalg.norm(w)
        
        ## Find 1st atom not colinear with 1 and 2
        nr = -1   # atomic number of noncolinear atom
        mindot = 1.0
        
        for i in range(2, natom):
            dn = dot(a_tmp_xyz[i], w_u)/linalg.norm(a_tmp_xyz[i])
            if fabs(dn) < mindot:
                mindot = dn
                nr = i
        #print dn 
        
        if nr == -1:
            myrot = tmp1rot
        elif not nr >= 2 and nr < natom:
            print "We're fucked"
        
        ### Determine 2nd rotation matrix
        
        if nr >= 2 and nr < natom:
           u = Txyz_gau[nr]
           v = a_tmp_xyz[nr]
           upw = dot(u, w_u)
           vpw = dot(v, w_u)
        
        ### Subract off projection onto w axis to make
        #  vectors orthogonal to rotation axis and normalize
        
           for i in range(3):
               u[i] = u[i] - upw*w_u[i]
               v[i] = v[i] - vpw*w_u[i]
           un = linalg.norm(u)
           vn = linalg.norm(v)
           for i in range(3):
               u[i] = u[i] / un
               v[i] = v[i] / vn
        
        ## Cross u and v for sin(th)
           tmp_vec = cross(u, v)
           sinth = linalg.norm(tmp_vec)
        
        ## Dot u and v to get cos(th)
           costh = dot(u, v)
        
        ## Find 2nd rotation matrix
        
           tmp2rot = RodRot(costh,sinth, w_u)
        ## Test back rotation for atom nr
           ext2 = zeros(3)
           ext2 = tmp_xyz[nr]
           M_ext2 = mat(ext2.copy())
           M_tmp2rot = mat(tmp2rot.copy())
           tmp2 = M_ext2 * M_tmp2rot
           a_tmp2 = numpy.asarray(tmp2)
        ## Check diff again
           mdiff = 0.0
           diff=zeros(3)
           for i in range(3):
               diff[i] = a_tmp2[0][i] - Txyz_gau[nr][i]
               if fabs(diff[i]) > mdiff:
                   mdiff = fabs(diff[i])
        
           #print mdiff
        
        ## If back rotation does not match, do reverse rotation
           if not mdiff < small:
               tmp2rot = RodRot(costh,sinth,w_r)
               M_tmp2rot = mat(tmp2rot.copy())
               tmp2 = M_ext2 * M_tmp2rot
               a_tmp2 = numpy.asarray(tmp2)
               mdiff = 0.0
               diff = zeros(3)
               for i in range(3):
                   diff[i] = a_tmp2[0][i] - Txyz_gau[nr][i]
                   if fabs(diff[i]) > mdiff:
                       mdiff = abs(diff[i])
         
               rotlogger.info("\nMax diff after 2nd rotation = %3.11f\n"
               %mdiff)
        
           myrot = M_tmp1rot * M_tmp2rot
        
        
        ## Rotate and translate back to orig geom and check
        
        Txyz = M_Txyz_ext * myrot
        a_Txyz = numpy.asarray(Txyz)
        xyz = zeros(natom*3)
        xyz.shape = (natom,3)
        mdiff=0.0
        for i in range(natom):
            for j in range(3):
                xyz[i][j] = a_Txyz[i][j] + xyz_gau[0][j]
                diff[j]=abs(xyz[i][j]-xyz_gau[i][j])
                if diff[j]>mdiff: mdiff=diff[j]

        rotlogger.info('\nAfter 3rd rotation: max diff = %3.11f\n' %mdiff)
# Begin 7/31 edit
        
        M_grad = mat(grad_ext.copy())
        RotGrad = M_grad * myrot
        rotlogger.info('\nRotation matrix:%s\n' %myrot)
        RotArray = numpy.asarray(RotGrad)
    
    else:
        RotArray = grad_ext
    if dip==True:
        return numpy.asmatrix(numpy.identity(3))
    return RotArray

 
def RodRot(cth,sth,vec):
    Rmat=numpy.zeros(9)
    Rmat.shape=(3,3)
    Rmat[0][0]=cth+vec[0]*vec[0]*(1-cth)
    Rmat[0][1] = -1*vec[2]*sth + vec[0]*vec[1]*(1-cth)
    Rmat[0][2] = vec[1]*sth + vec[0]*vec[2]*(1-cth)
    Rmat[1][0] = vec[2]*sth + vec[0] * vec[1]*(1-cth)
    Rmat[1][1] = cth + vec[1]*vec[1]*(1-cth)
    Rmat[1][2] = -1* vec[0]*sth + vec[1]*vec[2]*(1-cth)
    Rmat[2][0] = -1*vec[1]*sth + vec[0]*vec[2]*(1-cth)
    Rmat[2][1] = vec[0]*sth + vec[1]*vec[2]*(1-cth)
    Rmat[2][2] = cth + vec[2] *vec[2] * (1-cth)
    return Rmat
