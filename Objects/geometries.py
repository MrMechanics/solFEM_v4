#
#	geometries.py
#  -------------
#
#	This is the geometry module. It provides objects
#	needed to describe and manipulate 3D-geometry
#	for graphics. Including points, vectors,
#	lines, rotations and transformations. It is a modified
#	version built on code that was not written by me. 
#   I don't remember where I found the original code many
#   many years ago.
#


import math





class Point3D(object):
    def __init__(self,x=0.,y=0.,z=0.):
        self.coordinates = [x,y,z]
    def x(self):
        return self.coordinates[0]
    def y(self):
        return self.coordinates[1]
    def z(self):
        return self.coordinates[2]
    def __repr__(self):
        return "Point3D("+str(self.x())+","+str(self.y())+","+str(self.z())+")"
    def __str__(self):
        return "P("+str(self.x())+","+str(self.y())+","+str(self.z())+")"
    def get(self):
        return self.coordinates
    def returnCopy(self):
        return Point3D( self.x(), self.y(), self.z() )
    def asVector3D(self):
        return Vector3D( self.x(), self.y(), self.z() )
    def distance(self,other):
        return (other-self).length()
    def average(self,other):
        return Point3D( (self.x()+other.x())*0.5, (self.y()+other.y())*0.5, (self.z()+other.z())*0.5 )
    def __add__(self,other):
        return Point3D( self.x()+other.x(), self.y()+other.y(), self.z()+other.z() )
    def __sub__(self,other):
        if isinstance(other,Vector3D):
            return Point3D( self.x()-other.x(), self.y()-other.y(), self.z()-other.z() )
        return Vector3D( self.x()-other.x(), self.y()-other.y(), self.z()-other.z() )
    def __eq__(self,other):
        return self.x()==other.x() and self.y()==other.y() and self.z()==other.z()
    def __ne__(self,other):
        return not (self==other)





class Vector3D(object):
    def __init__(self,x=0.,y=0.,z=0.):
        self.coordinates = [x,y,z]
    def x(self):
        return self.coordinates[0]
    def y(self):
        return self.coordinates[1]
    def z(self):
        return self.coordinates[2]
    def __repr__(self):
        return "Vector3D("+str(self.x())+","+str(self.y())+","+str(self.z())+")"
    def __str__(self):
        return "V("+str(self.x())+","+str(self.y())+","+str(self.z())+")"
    def get(self):
        return self.coordinates
    def returnCopy(self):
        return Vector3D( self.x(), self.y(), self.z() )
    def asPoint3D(self):
        return Point3D( self.x(), self.y(), self.z() )
    def lengthSquared(self):
        return self.x()*self.x()+self.y()*self.y()+self.z()*self.z()
    def length(self):
        return math.sqrt( self.lengthSquared() )
    def normalized(self):
        l = self.length()
        if ( l > 0 ):
            return Vector3D( self.x()/l, self.y()/l, self.z()/l )
        return self.returnCopy()
    def __neg__(self):
        return Vector3D( -self.x(), -self.y(), -self.z() )
    def __add__(self,other):
        if isinstance(other,Point3D):
            return Point3D( self.x()+other.x(), self.y()+other.y(), self.z()+other.z() )
        return Vector3D( self.x()+other.x(), self.y()+other.y(), self.z()+other.z() )
    def __sub__(self,other):
        return Vector3D( self.x()-other.x(), self.y()-other.y(), self.z()-other.z() )
    def __mul__(self,other):
        if isinstance(other,Vector3D):
           # dot product
           return self.x()*other.x() + self.y()*other.y() + self.z()*other.z()
        # scalar product
        return Vector3D( self.x()*other, self.y()*other, self.z()*other )
    def __rmul__(self,other):
        return self*other
    def __div__(self,other):
        return Vector3D( self.x()/other, self.y()/other, self.z()/other )
    def __xor__(self,other):   # cross product
        return Vector3D(
            self.y()*other.z() - self.z()*other.y(),
            self.z()*other.x() - self.x()*other.z(),
            self.x()*other.y() - self.y()*other.x() )
    def __eq__(self,other):
        return self.x()==other.x() and self.y()==other.y() and self.z()==other.z()
    def __ne__(self,other):
        return not (self==other)
    def angle(self,other):
        if isinstance(other,Vector3D):
            v1_u = self.normalized()
            v2_u = other.normalized()
            angle = math.acos(v1_u*v2_u)
            if math.isnan(angle):
                if v1_u == v2_u:
                    return 0.0
                else:
                    return math.pi
            return angle





class Quaternion(object):
    '''
A quaternion object, with the option to add, subtract
and multiply with other quaternions.
'''
    def __init__(self, array = [0.0,0.0,0.0,0.0]):
        '''
    Initialize the quaternion object with the r, i, j and
    k values specified.
    '''
        self.r = array[0]
        self.i = array[1]
        self.j = array[2]
        self.k = array[3]
    def printQuat(self):
        '''
    Prints out the quaternion in standard format (and with
    angles?).
    '''
        print(str("%.3f" % self.r) + " " + \
              str("%.3f" % self.i) + "i " + \
              str("%.3f" % self.j) + "j " + \
              str("%.3f" % self.k) + "k")
    def add(self, quat):
        '''
    Add to another quaternion and return the result.
    '''
        return Quaternion([self.r + quat.r,
                           self.i + quat.i,
                           self.j + quat.j,
                           self.k + quat.k])
    def sub(self, quat):
        '''
    Subtract from another quaternion and return the result.
    '''
        return Quaternion([self.r - quat.r,
                           self.i - quat.i,
                           self.j - quat.j,
                           self.k - quat.k])
    def scale(self, scalar):
        '''
    Multiply quaternion with scalar and return the result.
    '''
        return Quaternion([self.r * scalar,
                           self.i * scalar,
                           self.j * scalar,
                           self.k * scalar])
    def multi(self, quat):
        '''
    Multiply with another quaternion and return the result.
    '''
        return Quaternion([self.r*quat.r - self.i*quat.i - self.j*quat.j - self.k*quat.k,
                           self.r*quat.i + self.i*quat.r + self.j*quat.k - self.k*quat.j,
                           self.r*quat.j - self.i*quat.k + self.j*quat.r + self.k*quat.i,
                           self.r*quat.k + self.i*quat.j - self.j*quat.i + self.k*quat.r])
    def conj(self):
        '''
    Gives the conjugate of the quaternion.
    '''
        return Quaternion([ self.r, -self.i, -self.j, -self.k])
    def inv(self):
        '''
    Invert quaternion.
    '''
        if (self.r == self.i == self.j == self.k == 0):
            print("Zero Quaternion has no inverse")
        else:
            scl = 1.0 / (self.r**2 + self.i**2 + self.j**2 + self.k**2)
            return Quaternion([ self.r * scl,
                               -self.i * scl,
                               -self.j * scl,
                               -self.k * scl])
    def norm(self):
        '''
    Gives the norm of the quaternion.
    '''
        return math.sqrt(self.r**2 + self.i**2 + self.j**2 + self.k**2)
    def unit(self):
        '''
    Gives the unit quaternion (versor).
    '''
        if (self.r == self.i == self.j == self.k == 0):
            print("Zero Quaternion has no versor")
        else:
            norm = self.norm()
            return Quaternion([self.r / norm,
                               self.i / norm,
                               self.j / norm,
                               self.k / norm])
    def recip(self):
        '''
    Gives the reciprocal of quaternion.
    '''
        if (self.r == self.i == self.j == self.k == 0):
            print("Zero Quaternion has no reciprocal")
        else:
            return self.conj().scale(1.0/(self.norm()**2))
    def vectorToQuat(self, v):
        '''
    Converts quaternion to a vector quaternion with the 
    same values as vector v.	
    '''
        self.r = 0.
        self.i = v.x()
        self.j = v.y()
        self.k = v.z()
    def axisAngleToQuat(self, rAx, rAng):
        '''
    Converts quaternion to a rotation quaternion with 
    axis rAx and angle rAng (in radians).
    '''
        rAx = rAx.normalized()
        theta = rAng/2.

        self.r = math.cos(theta)
        self.i = rAx[0]*math.sin(theta)
        self.j = rAx[1]*math.sin(theta)
        self.k = rAx[2]*math.sin(theta)
    def quatToAxisAngle(self):
        '''
    Gives the Axis and Angle of rotation quaternion.
    '''
        rAx = [self.i, self.j, self.k]
        rAng = math.arccos(self.r)*2
        return rAx, rAng
    def rotatePointAboutAxis(point,axis0,axis1,angle):
        '''
    First moves point, axis0 and axis1 together so that 
    axis0 is at the origin. Then creates the two quaternions 
    and rotates the point about the arbitrary axis defined
    by axis0 and axis1. Finally moves the point back so that 
    axis0 is at its original position again. 
    '''
        axis1_0 = axis1 - axis0
        point_0 = point - axis0
    
        qPnt0 = Quaternion()
        qPnt0.vectorToQuat(point_0)

        qRot = Quaternion()
        qRot.axisAngleToQuat(axis1_0.normalized(),angle)
        qPnt1 = qRot.multi(qPnt0.multi(qRot.conj()))

        return [qPnt1.i+axis0.x(), qPnt1.j+axis0.y(), qPnt1.k+axis0.z()]




class Line(object):
    def __init__(self,number,points):
        self.number = number
        self.points = []
        self.newPoints(points)
    def newPoints(self,points):
        for p in points:
            self.points.append(Point3D(p.x(),p.y(),p.z()))
    def length(self):
        pass
    
    
    


class Arc(Line):
    def __init__(self,number,points,radius,center,axis):
        super(Arc,self).__init__(number,points)
        self.radius = radius
        self.center = Point3D(center.x(),center.y(),center.z())
        self.axis = CoordSys3D(self.center,axis.x_vec,axis.y_vec)
    def length(self):
        pass
    
    
    


class Spline(Line):
    def __init__(self,number,points):
        super(Spline,self).__init__(number,points)
    def length(self):
        pass





class Edge(object):
    def __init__(self,lines):
        self.lines = lines

    



class Face(object):
    def __init__(self,edges):
        self.edges = edges
    def normal(self):
        pass

    



class CoordSys3D(object):
    def __init__(self,origin,x_vec,y_vec):
        self.origin = Point3D(origin.x(),origin.y(),origin.z())
        self.x_vec = x_vec.normalized()
        self.y_vec = y_vec.normalized()
        self.z_vec = self.x_vec.__xor__(self.y_vec)
        self.z_vec = self.z_vec.normalized()
        self.y_vec = self.z_vec.__xor__(self.x_vec)
    def mapFrom(self,other,point):
        # Map point from coordinates in other csys to coordinates in this csys
        if isinstance(other,CoordSys3D) and isinstance(point,Point3D):
            MT = Matrix4x4.translation( other.origin.asVector3D() - self.origin.asVector3D() )
            MR_other = Matrix4x4()
            MR_other.m = [ other.x_vec.x(), other.y_vec.x(), other.z_vec.x(), 0.0, 
                           other.x_vec.y(), other.y_vec.y(), other.z_vec.y(), 0.0,
                           other.x_vec.z(), other.y_vec.z(), other.z_vec.z(), 0.0,
                                       0.0,             0.0,             0.0, 1.0 ]
            MR_self = Matrix4x4()
            MR_self.m = [ self.x_vec.x(), self.y_vec.x(), self.z_vec.x(), 0.0, 
                          self.x_vec.y(), self.y_vec.y(), self.z_vec.y(), 0.0,
                          self.x_vec.z(), self.y_vec.z(), self.z_vec.z(), 0.0,
                                     0.0,            0.0,            0.0, 1.0 ]
            MR = MR_other*MR_self
            newCoords = MR*MT*point
            return newCoords




    
    

class Matrix4x4(object):
    def __init__(self):
        self.setToIdentity()
    def __str__(self):
        return str(self.m[0:4]) + "\n" + str(self.m[4:8]) + "\n" + str(self.m[8:12]) + "\n" + str(self.m[12:16])
    def get(self):
        return self.m
    def returnCopy(self):
        M = Matrix4x4()
        M.m = list(self.m)  # copy the list
        return M
    def setToIdentity(self):
        self.m = [ 1.0, 0.0, 0.0, 0.0,
                   0.0, 1.0, 0.0, 0.0,
                   0.0, 0.0, 1.0, 0.0,
                   0.0, 0.0, 0.0, 1.0 ]

    @staticmethod
    def translation( vector3D ):
        M = Matrix4x4()
        M.m[ 0] = 1.0;   M.m[ 4] = 0.0;   M.m[ 8] = 0.0;   M.m[12] = vector3D.x();
        M.m[ 1] = 0.0;   M.m[ 5] = 1.0;   M.m[ 9] = 0.0;   M.m[13] = vector3D.y();
        M.m[ 2] = 0.0;   M.m[ 6] = 0.0;   M.m[10] = 1.0;   M.m[14] = vector3D.z();
        M.m[ 3] = 0.0;   M.m[ 7] = 0.0;   M.m[11] = 0.0;   M.m[15] = 1.0;
        return M

    @staticmethod
    def rotationAroundOrigin( angleInRadians, axisVector ):
        # Note: assumes axisVector is normalized
        c = math.cos( angleInRadians )
        s = math.sin( angleInRadians )
        one_minus_c = 1-c
        M = Matrix4x4()
        M.m[ 0] = c + one_minus_c * axisVector.x()*axisVector.x()
        M.m[ 5] = c + one_minus_c * axisVector.y()*axisVector.y()
        M.m[10] = c + one_minus_c * axisVector.z()*axisVector.z()
        M.m[ 1] = M.m[ 4] = one_minus_c * axisVector.x()*axisVector.y();
        M.m[ 2] = M.m[ 8] = one_minus_c * axisVector.x()*axisVector.z();
        M.m[ 6] = M.m[ 9] = one_minus_c * axisVector.y()*axisVector.z();
        xs = axisVector.x() * s
        ys = axisVector.y() * s
        zs = axisVector.z() * s
        M.m[ 1] += zs;  M.m[ 4] -= zs;
        M.m[ 2] -= ys;  M.m[ 8] += ys;
        M.m[ 6] += xs;  M.m[ 9] -= xs;

        M.m[12] = 0.0;
        M.m[13] = 0.0;
        M.m[14] = 0.0;
        M.m[ 3] = 0.0;   M.m[ 7] = 0.0;   M.m[11] = 0.0;   M.m[15] = 1.0;
        return M

    @staticmethod
    def rotation( angleInRadians, axisVector, originPoint ):
        v = originPoint.asVector3D()
        return Matrix4x4.translation(v) * Matrix4x4.rotationAroundOrigin(angleInRadians,axisVector) * Matrix4x4.translation(- v)

    @staticmethod
    def uniformScaleAroundOrigin(scaleFactor):
        M = Matrix4x4()
        M.m[ 0] = scaleFactor; M.m[ 4] = 0.0;         M.m[ 8] = 0.0;         M.m[12] = 0.0;
        M.m[ 1] = 0.0;         M.m[ 5] = scaleFactor; M.m[ 9] = 0.0;         M.m[13] = 0.0;
        M.m[ 2] = 0.0;         M.m[ 6] = 0.0;         M.m[10] = scaleFactor; M.m[14] = 0.0;
        M.m[ 3] = 0.0;         M.m[ 7] = 0.0;         M.m[11] = 0.0;         M.m[15] = 1.0;
        return M

    @staticmethod
    def uniformScale( scaleFactor, originPoint ):
        v = originPoint.asVector3D()
        return Matrix4x4.translation(v) * Matrix4x4.uniformScaleAroundOrigin(scaleFactor) * Matrix4x4.translation(- v)

    @staticmethod
    def lookAt( eyePoint, targetPoint, upVector, isInverted ):
        # step one: generate a rotation matrix

        z = (eyePoint-targetPoint).normalized()
        y = upVector
        x = y ^ z   # cross product
        y = z ^ x   # cross product

        # Cross product gives area of parallelogram, which is < 1 for
        # non-perpendicular unit-length vectors; so normalize x and y.
        x = x.normalized()
        y = y.normalized()

        M = Matrix4x4()

        if isInverted :
            # the rotation matrix
            M.m[ 0] = x.x();   M.m[ 4] = y.x();   M.m[ 8] = z.x();   M.m[12] = 0.0;
            M.m[ 1] = x.y();   M.m[ 5] = y.y();   M.m[ 9] = z.y();   M.m[13] = 0.0;
            M.m[ 2] = x.z();   M.m[ 6] = y.z();   M.m[10] = z.z();   M.m[14] = 0.0;
            M.m[ 3] = 0.0;     M.m[ 7] = 0.0;     M.m[11] = 0.0;     M.m[15] = 1.0;

            # step two: premultiply by a translation matrix
            return Matrix4x4.translation( eyePoint.asVector3D() ) * M
        else:
            # the rotation matrix
            M.m[ 0] = x.x();   M.m[ 4] = x.y();   M.m[ 8] = x.z();   M.m[12] = 0.0;
            M.m[ 1] = y.x();   M.m[ 5] = y.y();   M.m[ 9] = y.z();   M.m[13] = 0.0;
            M.m[ 2] = z.x();   M.m[ 6] = z.y();   M.m[10] = z.z();   M.m[14] = 0.0;
            M.m[ 3] = 0.0;     M.m[ 7] = 0.0;     M.m[11] = 0.0;     M.m[15] = 1.0;

            # step two: postmultiply by a translation matrix
            return M * Matrix4x4.translation( - eyePoint.asVector3D() )

    def __mul__(a,b):   # note: a is really self
        if isinstance(b,Matrix4x4):
            M = Matrix4x4()
            M.m[ 0] = a.m[ 0]*b.m[ 0] + a.m[ 4]*b.m[ 1] + a.m[ 8]*b.m[ 2] + a.m[12]*b.m[ 3];
            M.m[ 1] = a.m[ 1]*b.m[ 0] + a.m[ 5]*b.m[ 1] + a.m[ 9]*b.m[ 2] + a.m[13]*b.m[ 3];
            M.m[ 2] = a.m[ 2]*b.m[ 0] + a.m[ 6]*b.m[ 1] + a.m[10]*b.m[ 2] + a.m[14]*b.m[ 3];
            M.m[ 3] = a.m[ 3]*b.m[ 0] + a.m[ 7]*b.m[ 1] + a.m[11]*b.m[ 2] + a.m[15]*b.m[ 3];

            M.m[ 4] = a.m[ 0]*b.m[ 4] + a.m[ 4]*b.m[ 5] + a.m[ 8]*b.m[ 6] + a.m[12]*b.m[ 7];
            M.m[ 5] = a.m[ 1]*b.m[ 4] + a.m[ 5]*b.m[ 5] + a.m[ 9]*b.m[ 6] + a.m[13]*b.m[ 7];
            M.m[ 6] = a.m[ 2]*b.m[ 4] + a.m[ 6]*b.m[ 5] + a.m[10]*b.m[ 6] + a.m[14]*b.m[ 7];
            M.m[ 7] = a.m[ 3]*b.m[ 4] + a.m[ 7]*b.m[ 5] + a.m[11]*b.m[ 6] + a.m[15]*b.m[ 7];

            M.m[ 8] = a.m[ 0]*b.m[ 8] + a.m[ 4]*b.m[ 9] + a.m[ 8]*b.m[10] + a.m[12]*b.m[11];
            M.m[ 9] = a.m[ 1]*b.m[ 8] + a.m[ 5]*b.m[ 9] + a.m[ 9]*b.m[10] + a.m[13]*b.m[11];
            M.m[10] = a.m[ 2]*b.m[ 8] + a.m[ 6]*b.m[ 9] + a.m[10]*b.m[10] + a.m[14]*b.m[11];
            M.m[11] = a.m[ 3]*b.m[ 8] + a.m[ 7]*b.m[ 9] + a.m[11]*b.m[10] + a.m[15]*b.m[11];

            M.m[12] = a.m[ 0]*b.m[12] + a.m[ 4]*b.m[13] + a.m[ 8]*b.m[14] + a.m[12]*b.m[15];
            M.m[13] = a.m[ 1]*b.m[12] + a.m[ 5]*b.m[13] + a.m[ 9]*b.m[14] + a.m[13]*b.m[15];
            M.m[14] = a.m[ 2]*b.m[12] + a.m[ 6]*b.m[13] + a.m[10]*b.m[14] + a.m[14]*b.m[15];
            M.m[15] = a.m[ 3]*b.m[12] + a.m[ 7]*b.m[13] + a.m[11]*b.m[14] + a.m[15]*b.m[15];

            return M
        elif isinstance(b,Vector3D):
            # We treat the vector as if its (homogeneous) 4th component were zero.
            return Vector3D(
                a.m[ 0]*b.x() + a.m[ 4]*b.y() + a.m[ 8]*b.z(), # + a.m[12]*b.w(),
                a.m[ 1]*b.x() + a.m[ 5]*b.y() + a.m[ 9]*b.z(), # + a.m[13]*b.w(),
                a.m[ 2]*b.x() + a.m[ 6]*b.y() + a.m[10]*b.z()  # + a.m[14]*b.w(),
              # a.m[ 3]*b.x() + a.m[ 7]*b.y() + a.m[11]*b.z()    + a.m[15]*b.w()
                )
        elif isinstance(b,Point3D):
            # We treat the point as if its (homogeneous) 4th component were one.
            return Point3D(
                a.m[ 0]*b.x() + a.m[ 4]*b.y() + a.m[ 8]*b.z() + a.m[12],
                a.m[ 1]*b.x() + a.m[ 5]*b.y() + a.m[ 9]*b.z() + a.m[13],
                a.m[ 2]*b.x() + a.m[ 6]*b.y() + a.m[10]*b.z() + a.m[14]
                )






