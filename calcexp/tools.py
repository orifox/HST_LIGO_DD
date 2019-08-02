#!/usr/bin/env python
import os, string, re, sys, math, types, pickle,six
from subprocess import Popen, PIPE, STDOUT
import datetime
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
    
try: 
    from configparser import SafeConfigParser
except:
    from ConfigParser import SafeConfigParser

from astropy.time import Time

class SMFieldSearch:
    def __init__(self, cmpin, ra, dec, path):
        self.cmpin = cmpin
        self.ra    = ra
        self.dec   = dec
        self.startpath = path

        self.photmatch = []
        self.cmpmatch  = []

        smfield, smamp = self.parseCmp(self.cmpin)
        if smfield and smamp:
            self.locateCmp(smfield, smamp)

        #print(self.ra, self.dec, self.cmpmatch)

        self.search()

        #print(self.photmatch)

    def search(self):
        for cmpMatch in self.cmpmatch:
            self.scanCmp(cmpMatch)

    def parseCmp(self, cmpfile):
        p = re.compile('(sm\d\d).\d{6,6}_\d{4,4}.\d{3,4}_(\d{1,2})_sub.cmp')
        m = p.search(cmpfile)
        if m:
            return m.group(1), m.group(2)
        else :
            return None, None

    def locateCmp(self, smfield, amp):
        matches = []
        for datedir in os.listdir(self.startpath):
            absdatedir = os.path.join(self.startpath, datedir)
            if os.path.isdir(absdatedir) and (string.find(datedir, '_') != -1):
                ampdir = os.path.join(absdatedir, str(amp))
                if os.path.isdir(ampdir):
                    for smfile in os.listdir(ampdir):
                        abssmfile = os.path.join(ampdir, smfile)
                        if smfile.startswith(smfield) and smfile.endswith('sub.cmp'):
                            self.cmpmatch.append(abssmfile)

    def scanCmp(self, cmpfile, pixrad=2):
        import pcfitsio

        if not os.path.isfile(cmpfile):
            return None

        xs, ys = list(map(float, sky2xy(cmpfile, self.ra, self.dec)))

        if xs == None:
            return None

        for line in open(cmpfile).readlines()[1:]:

            xc, yc = list(map(float, line.split()[0:2]))

            dist = math.sqrt( (xs-xc)**2 + (ys-yc)**2 )

            #print(self.ra, self.dec, cmpfile, xs, ys, xc, yc, dist)

            if dist < pixrad:
                cptr = pcfitsio.fits_open_file(cmpfile, 0)
                mjd  = pcfitsio.fits_read_key_dbl(cptr, 'MJD-OBS')[0]
                pcfitsio.fits_close_file(cptr)

                self.photmatch.append((mjd, line))

def sky2xy(fitsfile, ra, dec, coordsys='', opts=''):
    if not os.path.isfile(fitsfile):
        print(('ERROR sky2xy: could not find file ',fitsfile))
        return None, None
    so = os.popen('sky2xy %s %s %s %s %s' % (opts,fitsfile,str(ra),str(dec),coordsys))
    output = so.readlines()[0]
    so.close()

    return parse_sky2xy_output(output)

def sky2xy_list(fitsfile, radeclist, opts='-j',verbose=0):

    if not os.path.isfile(fitsfile):
        print(('error: fits file %s doesn\'t exist' % (fitsfile)))
        return 1,None

    cmd = 'sky2xy %s %s ' % (opts,fitsfile)
    for radec in radeclist: cmd += ' %s' % (radec)
    errorflag,output = executecommand(cmd,'',verbose=verbose)
    if len(output) != len(radeclist):
        print(('RA/DEC:',radeclist))
        print(('sky2xy output:',output))
        print('error: inconsistent number of output X/Y position to input RA/DEC positions!')
        return 2,None

    xydict = {}
    irange = list(range(len(output)))
    for i in irange:

        xydict[radeclist[i]] = {}
        xydict[radeclist[i]]['unparsed']  = output[i]

        xydict[radeclist[i]]['x'],xydict[radeclist[i]]['y'] = parse_sky2xy_output(output[i])
        if xydict[radeclist[i]]['x']==None or xydict[radeclist[i]]['y']==None:
            print(('ERROR! cannot parse sky2xy output:',output[i]))
            return 1, None

    del irange,output
    return 0,xydict

def sky2xy_file(fitsfile, radecfile, opts='',verbose=0):
    if not os.path.isfile(fitsfile):
        print(('error: fits file %s doesn\'t exist' % (fitsfile)))
        return 1,None

    cmd = 'sky2xy %s %s @%s' % (opts,fitsfile,radecfile)
    errorflag,output = executecommand(cmd,'')

    xylist = []
    irange = list(range(len(output)))
    for i in irange:
        x,y = parse_sky2xy_output(output[i])
        if x==None or y==None:
            print(('ERROR! cannot parse xy2sky output:',output[i]))
            return 1, None
        xylist.append((x,y))

    del irange,output
    return 0,xylist


def xy2sky(fitsfile, x, y, opts=''):
    if not os.path.isfile(fitsfile):
        return None, None

    so = os.popen('xy2sky ' + opts + ' ' + fitsfile + ' '+ str(x) + ' ' + str(y))
    output = so.readlines()[0]
    so.close()

    return parse_xy2sky_output(output)

def xy2sky_list(fitsfile, xylist, opts='-j'):

    if not os.path.isfile(fitsfile):
        print(('error: fits file %s doesn\'t exist' % (fitsfile)))
        return 1,None

    cmd = 'xy2sky %s %s ' % (opts,fitsfile)
    for xy in xylist: cmd += xy
    errorflag,output = executecommand(cmd,'')

    if len(output) != len(xylist):
        print(('XY:',xylist))
        print(('xy2sky output:',output))
        print('error: inconsistent number of output RA/DEC position to input X/Y positions!')
        return 2,None

    radecdict = {}
    irange = list(range(len(output)))
    for i in irange:

        radecdict[xylist[i]] = {}
        radecdict[xylist[i]]['unparsed']  = output[i]

        radecdict[xylist[i]]['ra'],radecdict[xylist[i]]['dec'] = parse_xy2sky_output(output[i])
        if radecdict[xylist[i]]['ra']==None or radecdict[xylist[i]]['dec']==None:
            print(('ERROR! cannot parse xy2sky output:',output[i]))
            return 1, None
    del irange,output
    return 0,radecdict

def parse_xy2sky_output(output):
    fields = output.split()
    if 'FK5' in fields:
        system = 'FK5'
    elif 'J2000' in fields:
        system = 'J2000'
    elif 'FK4' in fields:
        system = 'FK4'
    elif 'B1950' in fields:
        system = 'B1950'
    else:
        return None, None

    if fields.index(system) == 2:
        return fields[0], fields[1]
    elif fields.index(system) == 5:
        return fields[3], fields[4]
    else:
        return None, None

def parse_sky2xy_output(output):
    fields = output.split()

    if 'FK5' in fields:
        system = 'FK5'
    elif 'J2000' in fields:
        system = 'J2000'
    elif 'FK4' in fields:
        system = 'FK4'
    elif 'B1950' in fields:
        system = 'B1950'
    else:
        return None, None

    if '(off' in fields or 'off' in fields:
        return None, None

    if fields.index(system) == 2:
        return list(map(float, (fields[4], fields[5])))
    elif fields.index(system) == 5:
        return list(map(float, (fields[3], fields[4])))
    else:
        return None, None

def calcPA(ra1,dec1,ra2,dec2):
    deg2rad=0.0174532925199      # math.pi / 180.
    rad2deg=57.2957795131      # 180.0 / math.pi
    ra1=RaInDeg(ra1)
    ra2=RaInDeg(ra2)
    dec1=DecInDeg(dec1)
    dec2=DecInDeg(dec2)
    dec1rad=dec1*deg2rad
    dec2rad=dec2*deg2rad

    if ra1-ra2>180:
        ra2=ra2+360.0
    if ra1-ra2<-180:
        ra2=ra2-360.0

    dRa = ra2*deg2rad - ra1*deg2rad

    pa_deg = rad2deg*math.atan2(math.sin(dRa),math.cos(dec1rad)*math.sin(dec2rad)/math.cos(dec2rad)-math.sin(dec1rad)*math.cos(dRa));
    return pa_deg

# this uses skycoor from WCStools, but there is a bug in skycoor ...
def calcPAbad(ra1,dec1,ra2,dec2):
    ra1=RaInDeg(ra1)
    ra2=RaInDeg(ra2)
    dec1=DecInDeg(dec1)
    dec2=DecInDeg(dec2)
    # This is to get rid of a skycoor bug:
    # rest@000d60780e33(SM,v9.0)% skycoor -a 350.85 58.815 4.66116 61.87606
    #-66.791
    #rest@000d60780e33(SM,v9.0)% skycoor -a -9.15 58.815 4.66116 61.87606
    #66.791
    if ra1-ra2>180:
        ra2=ra2+360.0
    if ra1-ra2<-180:
        ra2=ra2-360.0

    #cmd = 'skycoor -a %f %f %f %f' % (ra1,dec1,ra2,dec2)
    cmd = 'skycoor -a %s %s %s %s' % (deg2sex(ra1,ra=True),deg2sex(dec1),deg2sex(ra2,ra=True),deg2sex(dec2))
    #print(cmd)
    errorflag,skycoorout = executecommand(cmd,'',verbose=0)
    if errorflag or len(skycoorout)<1:
        print(('ERROR: could not determine PA!',skycoorout))
        PA = None
    else:
        PA = float(skycoorout[0])
        if PA<0.0: PA+=360.0
    return PA

# this routine moves length_arcsec along the position angle PA_deg from RA, DEC
# PA:
# 0..90 degree: NE quadrant
# 90..180:      SE quadrant
# 0..-90:       NW quadrant
# -90..-180:    SW quadrant
# since skycoor has only accuracy to 1/1000 of degree, we can get the PA only to 3.6 arcsec
def moveRADEC(ra, dec, PA_deg, length_arcsec, precision_PA_arcsec = 3.6, precision_length_arcsec = 0.001, precision_arcsec=0.0001, verbose=0):
    deg2rad = 0.0174532925199  # math.pi / 180.0
    rad2deg = 57.2957795131    # 180.0 / math.pi
    arcsec2deg = 1.0/3600.0
    deg2arcsec = 3600.0

    PA_rad = PA_deg*deg2rad
    precision_PA_deg = precision_PA_arcsec * arcsec2deg
    precision_length_deg = precision_length_arcsec * arcsec2deg
    precision_deg = precision_arcsec * arcsec2deg
    length_deg = length_arcsec * arcsec2deg

    RAdeg = sex2deg(ra,ra=True)
    DECdeg = sex2deg(dec)
    if abs(DECdeg)>=90.0:
        print('ERROR: abs(Dec)>=90.0 are not allowed!')
        sys.exit(0)
    stopflag = 0
    cosDEC = math.cos(DECdeg*deg2rad)
    deltaRA = length_deg * math.sin(PA_rad)/cosDEC
    deltaDEC = length_deg * math.cos(PA_rad)
    iteration = 0
    if verbose>=2:
        print(('RA:%f DEC:%f PA:%.5f deg length:%.8f deg' % (RAdeg,DECdeg,PA_deg,length_deg)))
    while not stopflag:

        if verbose>=2:
            distance0 = skydist_degree(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)

        PAtest_deg = calcPA(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)
        if PAtest_deg-PA_deg>=90.0:PAtest_deg-=180.0
        if PAtest_deg-PA_deg<=-90.0:PAtest_deg+=180.0

        # First adjust the PA
        dDec = -length_deg*math.sin(PAtest_deg*deg2rad)*(PA_rad - PAtest_deg*deg2rad)
        dRa  = length_deg*math.cos(PAtest_deg*deg2rad)*(PA_rad - PAtest_deg*deg2rad)/cosDEC

        deltaDEC += dDec
        deltaRA  += dRa
        if abs(DECdeg+deltaDEC) >= 90:
            print('ERROR: the position is too close to the pole, not yet implemented!')
            return(2,None,None)

        if verbose>=2:
            PAtest1_deg = calcPA(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)
            if PAtest1_deg-PA_deg>=90.0:PAtest1_deg-=180.0
            if PAtest1_deg-PA_deg<=-90.0:PAtest1_deg+=180.0
            distance1 = skydist_degree(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)

        # Now adjust the length
        distance = skydist_degree(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)
        distance_ratio = length_deg/distance
        deltaDEC *= distance_ratio
        deltaRA *= distance_ratio
        distance = skydist_degree(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)

        if verbose>=2:
            PAtest2_deg = calcPA(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)
            if PAtest2_deg-PA_deg>=90.0:PAtest2_deg-=180.0
            if PAtest2_deg-PA_deg<=-90.0:PAtest2_deg+=180.0
            distance2 = skydist_degree(RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)

        iteration +=1
        goodflag = (abs(length_deg - distance) < precision_length_deg) and (abs(PAtest_deg - PA_deg)<precision_PA_deg) \
            or abs(dRa)+abs(dDec)<precision_deg

        stopflag = (iteration>=40) or goodflag

        if verbose==1:
            print(('iteration %d: %.8f %.8f -> %.8f %.8f   (%.8e %.8e)' % (iteration,RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC,dRa,dDec)))
        elif verbose >=2:
            print(('iteration %d: PA0:%.8f PA1:%.8f PA2:%.8f d0:%.8f d1:%.8f d2:%.8f deltaRa:%11.8f deltaDec:%11.8f (%.4e %.4e)' % (iteration,PAtest_deg,PAtest1_deg,PAtest2_deg,distance0,distance1,distance2,deltaRA,deltaDEC,dRa,dDec)))
            print(('(%.8e %.8e %.8e %.8e %.8e %.8e)' % (abs(length_deg - distance),precision_length_deg,abs(PAtest_deg - PA_deg),precision_PA_deg,abs(dRa)+abs(dDec),precision_deg)))


    #print('results:',RAdeg,DECdeg,RAdeg+deltaRA,DECdeg+deltaDEC)
    return(not goodflag,RAdeg+deltaRA,DECdeg+deltaDEC)


### Converts sexigesimal notation to decimal degrees or decimal hours (if option 'ra=True' invoked)
def sex2deg(sexigecimal, ra=False):
    print(sexigecimal)
    ### 2005/12/02 - AR: make sure it is in sexagesimal format!
    # is it a string? if not check if it is None
    if not (type(sexigecimal) is str): #types.StringType):
        if type(sexigecimal) == None:
            raise RuntimeError("ERROR: sex2deg cannot handle 'None'!")
        return sexigecimal
    # Does it have ':' or ' '? If not, it must be a float in string format, just convert it and return
    if re.search('\:|\s',sexigecimal) == None:
        return(float(sexigecimal))

    s1, s2, s3 = list(map(float, re.split('[DHMSdhms:\s]', sexigecimal.strip())[0:3]))
    # Get the sign
    if re.search('-', sexigecimal):
        sign = -1
    else:
        sign = 1

    deg = abs(s1) + s2 / 60. + s3 / 3600.

    deg *= sign

    if ra:
        deg *= 15.

    return deg

### Converts decimal degrees or hours (ra=True) to sexigesimal notation
###  [-+]DD:MM:SS.ss
### the format is fixed at two decimals of precision for the decimal seconds.
### No -/+ if ra=True
def deg2sex(degrees, ra=False, outputformatRA='%02d:%02d:%06.3f',outputformatDEC='%1s%02d:%02d:%05.2f'):
    if type(degrees) is str: # types.StringType:
        degrees=float(degrees)
    if ra:
        # a.k.a. indeg and outhours
        if degrees < 0.0:
            while degrees<0:degrees+=360.
        if degrees > 360.0:
            while degrees>360.0:degrees-=360.
        degrees /= 15.

    if degrees < 0:
        sign = '-'
    else:
        sign = '+'

    degrees = abs(degrees)

    d1  = (degrees - (degrees % 1))
    rem = (degrees % 1) * 60
    d2  = rem - (rem % 1)
    srem = (rem % 1) * 60
    d3 = srem

    if ra:
      return outputformatRA % (d1, d2, d3)
    else:
      return outputformatDEC % (sign, d1, d2, d3)

# Returns the passed in RA in decimal degrees
# input RA can be in 'HH:MM:SS.ss', 'HH MM SS.ss' or in decimal degrees
def RaInDeg(Ra):
    import types
    if type(Ra)==str: #types.StringType:
        if re.search('[DHMShms: ]',Ra.strip()):
            return(sex2deg(Ra,ra=True))
        return(float(Ra))
    else:
        return(Ra)

# Returns the passed in Dec in decimal degrees
# input Dec can be in 'DD:MM:SS.ss', 'DD MM SS.ss' or in decimal degrees
def DecInDeg(Dec):
    import types
    if type(Dec)==str:#types.StringType:
        if re.search('[DHMShms: ]',Dec.strip()):
            return(sex2deg(Dec,ra=False))
        return(float(Dec))
    else:
        return(Dec)

def executecommand(cmd,successword,errorlog=None,cmdlog=None,verbose=1):
    if verbose: print(('executing: ',cmd))

#    (cmd_in,cmd_out)=os.popen4(cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    cmd_in,cmd_out = p.stdin,p.stdout

    output = cmd_out.readlines()
    if successword=='':
        successflag = 1
    else:
        m = re.compile(successword)
        successflag = 0
        for line in output:
            if sys.version_info[0] >= 3:
                line = line.decode('utf-8')
            if m.search(line):
                successflag = 1
    errorflag = not successflag
    if errorflag:
        print(('error executing:',cmd))
        if errorlog != None:
            append2file(errorlog,['\n error executing: '+cmd+'\n'])
            append2file(errorlog,output)
        if cmdlog != None:
            append2file(cmdlog,['\n error executing:',cmd])
    else:
        if cmdlog != None:
            append2file(cmdlog,[cmd])
    return errorflag, output

def append2file(filename,lines,verbose=0):
    if type(lines) is str:#types.StringType:
        lines = [lines,]
    if os.path.isfile(filename):
        buff = open(filename, 'a')
    else:
        buff = open(filename, 'w')
    r=re.compile('\n$')
    for line in lines:
        if not r.search(line):
            buff.write(line+'\n')
            if verbose: print((line+'\n'))
        else:
            buff.write(line)
            if verbose: print(line)
    buff.close()

def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=1):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

def hex2int(val):
    if type(val) is str:#types.StringType:
        val = int(eval(val))
    return(val)

def rmfile(filename,raiseError=1,gzip=False):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    if gzip and os.path.lexists(filename+'.gz'):
        os.remove(filename+'.gz')
        if os.path.isfile(filename+'.gz'):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename+'.gz')
            else:
                return(2)
    return(0)

def rmfiles(filenames,raiseError=1,gzip=False):
    if not (type(filenames) is list):
        raise RuntimeError("List type expected as input to rmfiles")
    errorflag = 0
    for filename in filenames:
        errorflag |= rmfile(filename,raiseError=raiseError,gzip=gzip)
    return(errorflag)

# some little defs to compare lists and tuples
def unique(seq):
    if seq == None or seq == []: return []
    d = {}
    for x in seq:
        d[x] = 1
    return list(d.keys())

# some little defs to compare lists and tuples
def unique_keeporder(seq):
    if seq == None or seq == []: return []
    d = {}
    for x in seq:
        d[x] = 0
    list = []
    for x in seq:
        d[x] += 1
        if d[x]==1:
            list.append(x)

    return list

def multiple(seq):
    if seq == None: return []
    d = {}
    for x in seq:
        if x in d:
            d[x] += 1
        else:
            d[x] = 1
    m = []
    for x in list(d.keys()):
        if d[x]>1:
            m.append(x)
    return m

def AnotB(A,B):
    "returns elements that are in A, but not in B"
    if A == None: return []
    if not (type(A) in [list,tuple]): A = [A,]
    if B == None: return A
    if not (type(B) in [list,tuple]): B = [B,]
    c = {}
    d = {}
    for x in A:
        c[x] = 1
        d[x] = 1
    for x in B:
        if x in c:
            if x in d:
                del(d[x])
    del c
    return list(d.keys())

def AorB(A,B):
    "returns list of elements that are in A or in B"
    if A == None and B == None: return []
    if B != None:
        if not (type(B) in [list,tuple]): B = [B,]
    if A != None:
        if not (type(A) in [list,tuple]): A = [A,]
    if A == None: return B
    if B == None: return A

    d = {}
    for x in A:
        d[x] = 1
    for x in B:
        d[x] = 1
    return list(d.keys())

#def AandB(A,B):
#    "returns list of elements that are in A and in B"
#    if A == None or B == None: return []
#    if not (type(A) in [types.ListType,types.TupleType]): A = [A,]
#    if not (type(B) in [types.ListType,types.TupleType]): B = [B,]
#    c = {}
#    d = {}
#    for x in A:
##        c[x] = 1
#    for x in B:
#        if c.has_key(x):
#            d[x] = 1
#    del c
#    return d.keys()

def AandB(A,B):
    "returns list of elements that are in A and in B, keeps the order from A"
    if A == None or B == None: return []
    if not (type(A) in [list,tuple]): A = [A,]
    if not (type(B) in [list,tuple]): B = [B,]
    c = {}
    for x in B:
        c[x] = 0
    ablist = []
    for x in A:
        if x in c:
            ablist.append(x)
        c[x] = 1
    del c
    return ablist

def not_AandB(A,B):
    "returns elements that are only in A or only in B"
    if A == None or B == None: return self.AorB(A,B)
    if not (type(A) in [list,tuple]): A = [A,]
    if not (type(B) in [list,tuple]): B = [B,]
    d = {}
    for x in B:
        d[x] = 1
    for x in A:
        d[x] = 1
    crosssection = AandB(A,B)
    for x in crosssection:
        if x in d:
            del(d[x])
    return list(d.keys())

#approximate sky dist (good for small angular distances not too close to pole)
def approxskydist_degree(ra1,dec1,ra2,dec2):
    deg2rad=0.0174532925199  # math.pi / 180.
    #convert all values to radians
    ra1 =RaInDeg(ra1)
    ra2 =RaInDeg(ra2)
    dec1=DecInDeg(dec1)
    dec2=DecInDeg(dec2)
    cosdec=min(abs(math.cos(dec1*deg2rad)),abs(math.cos(dec2*deg2rad)))
    return (math.sqrt(cosdec*cosdec*(ra2-ra1)*(ra2-ra1)+(dec2-dec1)*(dec2-dec1)))

# true angular distance between two position in the sky
#http://www2.sjsu.edu/faculty/watkins/sphere.htm
#cos(A) = sin(phi2)sin(phi1)+ cos(phi2)cos(phi1)cos(theta2 - theta1)
#longitude=theta=RA
#latitude=phi=DEC
#cos(A) = sin(DEC2)sin(DEC1)+ cos(DEC2)cos(DEC1)cos(RA2 - RA1)
#works also at RA=0 etc:
#skydist_degree(ra1,dec1,ra2,dec2) = skydist_degree(ra1+360,dec1+360,ra2-360,dec2-360) etc!
def skydist_degree(ra1,dec1,ra2,dec2):
    deg2rad=0.0174532925199  # math.pi / 180.
    #convert all values to radians
    ra1 =RaInDeg(ra1)*deg2rad
    ra2 =RaInDeg(ra2)*deg2rad
    dec1=DecInDeg(dec1)*deg2rad
    dec2=DecInDeg(dec2)*deg2rad
    if ((ra1==ra2) and (dec1==dec2)): return 0
    return math.acos(math.sin(dec2)*math.sin(dec1)+ math.cos(dec2)*math.cos(dec1)*math.cos(ra2 - ra1))*180.0/math.pi

def sphcor(decdeg):
    # spherical correction for RA distances as a function of dec
    # as in delta_radeg *= sphcor(decdeg)
    return math.cos(decdeg * math.pi / 180.)

# returns a dictionary d[file][fitskey] containing the values of the fitskeys,
# with file and fitskey from the filelist and keylist
# if a file doesn't exist, d[file]=None, or exception raised if erroriffilenotexist
# if a fitskey doesn't exist, then d[file][fitskey] == None, or exception raised if errorifkeynotexist
def getfitskeys(filelist,fitskeylist,erroriffilenotexist=True,errorifkeynotexist=False):
    if isinstance(filelist,six.string_types): filelist     = [filelist,]
    if isinstance(fitskeylist,six.string_types): fitskeylist  = [fitskeylist,]
    if not (type(fitskeylist)  is list):
        raise RuntimeError('Error: def getfitskeys expects a list of fits keys')
    if not (type(filelist)     is list):
        raise RuntimeError('Error: def getfitskeys expects a list of files')

    # GSN - 20100615 - replace -pau with -paub so that blanks in the fits keys don't cause a bork
    cmd = 'gethead -paub '
    filelist2use = []
    filelist2skip = []
    for file in filelist:
        if os.path.isfile(file):
            cmd = cmd + ' %s' % file
            filelist2use.append(file)
        else:
            if erroriffilenotexist:
                raise RuntimeError('ERROR: file %s does not exist!' % (file))
            filelist2skip.append(file)

    for fitskey in fitskeylist:
        cmd = cmd + ' %s' % fitskey
    #(gethead_in,gethead_out)=os.popen4(cmd)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    gethead_in,gethead_out = p.stdin,p.stdout
    lines = gethead_out.readlines()

    if len(lines)!=len(filelist2use):
        for line in lines: print(line)
        raise RuntimeError('ERROR: something went wrong with '+cmd+'! output has different # of lines than there are files!')

    r = re.compile('[A-Z]+:')
    fitsvaldict = {}
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
        entries = re.split('\s+',lines[i].decode('utf-8'))
        # error checking
        if not (entries[0] == filelist2use[i]):
            raise RuntimeError('ERROR: something went wrong with %s! output has file %s, which is not equal to %s' % (cmd,entries[0],filelist2use[i]))
        if len(entries)-1 != len(fitskeylist):
            print(('parsed line:',entries))
            raise RuntimeError('ERROR: something went wrong with %s! file %s returned different number of values (%d) then fits keys (%d): %s' % (cmd,filelist2use[i],len(entries)-1,len(fitskeylist),lines[i]))
        if r.match(entries[1]):
            raise RuntimeError('ERROR: something went wrong with %s! There seems to be a gethead error message (%s) for file %s:' % (cmd,entries[1],filelist2use[i]))

        # fill the dictionary
        fitsvaldict[filelist2use[i]]={}
        for f in range(len(fitskeylist)):
            if entries[f+1] == '___':
                if errorifkeynotexist: # raise exception if wanted
                    raise RuntimeError('ERROR: key %s does not exist in file %s!' % (fitskeylist[f],filelist2use[i]))
                entries[f+1] = None
            fitsvaldict[filelist2use[i]][fitskeylist[f]] = entries[f+1]

    # set the files that don't exist to None
    for file in filelist2skip: fitsvaldict[file] = None

    return(fitsvaldict)

def check4default(val,default):
    if val=="default":
        returnval=default
    else:
        returnval=val
    return(returnval)

# returns c, so that f_ref = f *c, where f is the flux for a zeropoint
# of ZPT, and f_ref is the corresponding flux at a zeropoint of ZPTref
def c_f2fref(ZPTref,ZPT):
    c = 10.0**(-0.4*(ZPT - ZPTref))
    return(c)

def mag2flux(mag, dmag, offset=0):
    flux  = 10.**(-0.4 * (mag-offset))
    dflux = 0.4 * math.log(10.0) * dmag * flux
    return flux, dflux

def flux2mag(flux, dflux):
    mag  = -2.5 * log10(flux)
    dmag =  2.5 / log(10.0) * dflux / flux
    return mag, dmag

def getkeys(hdr, comp):
    keys = []
    cstr = re.compile(comp)
    keys = hdr.getkeys(comp)
    return keys

def calcMJD(datetimestructure):
    # 1858-11-17 00:00:00.00 = MJD 0
    date0   = datetime.datetime(1858, 11, 17, 0, 0, 0, 0)
    delta = datetimestructure - date0
    MJD =  delta.days + (delta.seconds + delta.microseconds / 1.e6) / 86400.
    return(MJD)

def returnMJD():
    import datetime
    now    = datetime.datetime.now()

    # 1858-11-17 00:00:00.00 = MJD 0
    MJD0   = datetime.datetime(1858, 11, 17, 0, 0, 0, 0)

    delta = now - MJD0
    return delta.days + (delta.seconds + delta.microseconds / 1.e6) / 86400.

def pickledfilename(filename):
    if test4pickled(filename):
        return(filename)
    else:
        return(filename+'.pickled')

def unpickledfilename(filename):
    if test4pickled(filename):
        filename = re.sub('\.pickled$','',filename)
        return(filename)
    else:
        return(filename)

def test4pickled(filename):
    if re.search('\.pickled$',filename):
        return(1)
    else:
        return(0)

def pickleobject(object, outfile):
    f = open(outfile, 'w')
    p = pickle.Pickler(f)
    p.dump(object)
    f.close()

def unpickleobject(filename):
    f = open(filename, 'r')
    p = pickle.Unpickler(f)
    object = p.load()
    f.close()
    return object

def isfile_or_gz(filename,checkfilesize=False,checkdouble=False):
    isfileFlag = os.path.isfile(filename)
    if os.path.isfile(filename):
        if checkfilesize and os.path.getsize(filename)==0:
            return(False)
        if checkdouble and os.path.isfile(filename+'.gz'):
            return(False)
        return(True)
    if os.path.isfile(filename+'.gz'):
        if checkfilesize and os.path.getsize(filename+'.gz')==0:
            return(False)
        return(True)
    return(False)

def SetFitsKeywords(filename,args,options='',verbose=False):
    cmd = 'sethead %s %s %s 2>&1' % (filename,options,args)
    if verbose: print(cmd)
    so = os.popen(cmd)
    output = so.readlines()
    so.close()
    if len(output)>0:
        print(('ERROR: something went wrong setting the fits keywords with command %s:' % cmd))
        for line in output:
            print(line)
        return 1
    return 0

def GetFitsKeywords(filename,keywordlist,options='',verbose=False):
    if (not os.path.isfile(filename)):
        return 10,[]

    cmd = 'gethead -g %s %s' % (options,filename)
    for keyword in keywordlist:
        cmd += ' '+keyword
    cmd += ' 2>&1'
    if verbose: print(cmd)
    so = os.popen(cmd)
    output = so.readlines()
    so.close()
    if len(output)!=len(keywordlist):
        print(('ERROR: something went wrong getting the fits keywords with command %s, inconsistent number of output lines %d!=%d' % (cmd,len(output),len(keywordlist))))
        for line in output:
            print(line)
        return 11,[]
    errorflag = 0
    keyvals = []
    for i in range(len(output)):
        m = re.match('(%s)\s*=\s*(.*)' % keywordlist[i],output[i])
        if m and len(m.groups())==2:
            keytest = m.groups()[0]
            if keytest != keywordlist[i]:
                print(('ERROR: inconsistent keyword returned from gethead %s!=%s' % (keywordlist[i],keytest)))
                return 12,[]
            keyvals.append(m.groups()[1])
        else:
            if re.match('%s\s+not\s+found' % keywordlist[i],output[i]):
                keyvals.append(None)
                errorflag = 1
            else:
                print(('ERROR: something is wrong with keyword %s: %s' % (keywordlist[i],output[i])))
                return 13,[]
    return errorflag,keyvals

# writeFits*: from actions.py

def writeFitsStamp(lc, lckey, fname, outfile, boxsizePix=200):
    fitsfile = lc.getentry(lckey, fname)
    if fitsfile == None:
        return 1

    ra       = lc.getentry(lckey, 'Ra')
    dec      = lc.getentry(lckey, 'Dec')
    chopstr  = ','.join((ra, dec, str(boxsizePix)+'p'))
    cmdStr   = ' '.join(('fitscopy -w', chopstr, fitsfile, outfile))
    print(cmdStr)
    cmd      = command.Command(cmdStr)
    cmd.execute()

    if not os.path.isfile(outfile):
        print(('Problem saving', outfile))
        return 1

    return 0


def writeFitsStamps(lc, lckey, saveImAs=None, saveTmplAs=None, saveSubAs=None, boxsizePix=200):
    if saveImAs != None:
        writeFitsStamp(lc, lckey, 'longinimfitsfile', saveImAs,   boxsizePix)
    if saveTmplAs != None:
        writeFitsStamp(lc, lckey, 'longtmplfitsfile', saveTmplAs, boxsizePix)
    if saveSubAs != None:
        writeFitsStamp(lc, lckey, 'longsubfitsfile',  saveSubAs,  boxsizePix)

# Armin: HACK: temporary, I have to pass directly the fits filename
def writeFitsStampTemporary(fitsfile, ra, dec, outfile, boxsizePix=200):
    chopstr  = ','.join((ra, dec, str(boxsizePix)+'p'))
    cmdStr   = ' '.join(('fitscopy -w', chopstr, fitsfile, outfile))
    print(cmdStr)
    (cmd_in,cmd_out) = os.popen4(cmdStr)
    for line in cmd_out.readlines():
        print(('Potential problem saving', outfile,':',line))
        pass

    if not os.path.isfile(outfile):
        print(('Problem saving', outfile))
        return 1

    return 0

def writeAllImages(outname, fitsfiles, badimflag, ra, dec, size=150, nX=8, xhair=1,MJDbase=0.0,xflipflag=0):
    import imDisp
    from PIL import Image, ImageDraw, ImageFont
    #import numarray as num
    import numpy as num
    #import mx.DateTime

    print(('VVVVVVVVVVVVVVV',xflipflag))

    phot2filters= {0x00:'VR',0x02:'B',0x03:'V',0x04:'R',0x05:'I',
                   0x12:'u',0x13:'g',0x14:'r',0x15:'i',0x16:'z',0x17:'y'}
                   

    font  = 'helvB08.pil'
    #dfont = ImageFont.load_path(font)
    dfont = None

    nX    = min(nX, len(fitsfiles))

    #nY      = int(len(fitsfiles) / nX) + 1
    # Armin: above made two rows if nX=3 and len(fitsfiles)=3
    nY      = int(len(fitsfiles) / nX)
    if len(fitsfiles) / float(nX) > nY:
        nY += 1

    buff1   = 2
    sizeX   = nX * size + (nX - 1) * buff1
    sizeY   = nY * size + (nY - 1) * buff1
    #comp    = Image.new('L', (sizeX, sizeY), 255)
    comp    = Image.new('RGB', (sizeX, sizeY), "white")

    X = 0
    Y = 0
    n = 0

    data = num.zeros((size, size))


    #print("FITSFILES: ", fitsfiles)
    for i in range(len(fitsfiles)):
        fitsfile = fitsfiles[i]
        if fitsfile == None:
            continue

        if not os.path.isfile(fitsfile):
            continue

        try:
            fptr = pyfits.open(fitsfile)
        except:
            print('ERROR: could not open fitsfile')
            continue

        ##
        # get header and scaling
        if fitsfile.endswith('diff.fits'):
            # these stats are not calculated anymore...
            skykey  = 'DMEAN00'
            sigkey  = 'DSIGE00'
        else:
            skykey = 'SKYADU'
            sigkey = 'SKYSIG'

#	try:
#            # D. Jones - adjusting for new pyfits version
#	    sky = float(skykey in fptr[0].header)
#	    sky = float(fptr[0].header.has_key(skykey))
#	except:
#	    print('WACKY ERROR.. file exists but key index odd...',fitsfile)
#	    continue

        if skykey in fptr[0].header:
#        if fptr[0].header.has_key(skykey):
            sky      = float(fptr[0].header[skykey])
            #print('Found sky in header key %s: %.1f' % (skykey,sky))
        else:
            print(('WARNING! could not find sky %s in header of %s' % (skykey,fitsfile)))
            sky      = 0.1

        if sigkey in fptr[0].header:
#        if fptr[0].header.has_key(sigkey):
            sigsky   = float(fptr[0].header[sigkey])
            #print('Found skysig in header key %s: %.1f' % (sigkey,sigsky))
        else:
            sigsky   = num.sqrt(sky)
            print(('WARNING! could not find skysig %s in header of %s' % (sigkey,fitsfile)))

        #print('Image %s: SKYADU:%.2f  SKYSIG:%.2f' % (fitsfile,sky,sigsky))

        if 'FILTER' in fptr[0].header:
            filter  = fptr[0].header['FILTER'][0:2]
        else:
            if 'PHOTCODE' in fptr[0].header:
                phot = hex2int(fptr[0].header['PHOTCODE']) & 0xff
                if phot in phot2filters:
                    filter  = phot2filters[phot]
                else:
                    filter='n/a'
        datestr = None
        if 'MJD-OBS' in fptr[0].header:
            dateobject=Time(fptr[0].header['MJD-OBS'],format='mjd', scale='utc')
            datestr=dateobject.iso  
        elif 'DATE-OBS' in fptr[0].header:
#        if fptr[0].header.has_key('DATE-OBS'):
            datestr = fptr[0].header['DATE-OBS']
#        elif 'DATE' in fptr[0].header:
#        elif fptr[0].header.has_key('DATE'):
#            datestr  = fptr[0].header['DATE']
        #try:
        #    date = mx.DateTime.DateTimeFrom(re.sub('[\'T]', ' ', datestr))
        #    datestr = date.date
        #    mjdstr  = date.mjd
        #except:
        #    # This is a hack: DateTimeFrom does not take 60.0 seconds
        #    temp = re.sub('60\.0$', '59.99999999', datestr)
        #    date = mx.DateTime.DateTimeFrom(re.sub('[\'T]', ' ', temp))
        #    datestr = date.date
        #    mjdstr  = date.mjd
        #label   = "%s : %8.2f : %s" % (date.date, date.mjd-MJDbase, filter)

        diffimflag=0
        if 'TMPLSUBD' in fptr[0].header:
#        if fptr[0].header.has_key('TMPLSUBD'):
            diffimflag=1

        if datestr != None:
            #label   = "%s : %s" % (datestr, filter)
#            label   = "%s :" % (datestr, filter)
#            flabel   = "%s" % (filter)
            label   = " : %s" % (datestr)
            if diffimflag:
                flabel   = "%s" % (filter)
            else:
                flabel   = "%s" % (filter)
        else:
            #label   = "DATE UNKOWN: %s" % (filter)
            label   = " : DATE UNKOWN"
        xfits   = fptr[0].header['NAXIS1']
        yfits   = fptr[0].header['NAXIS2']

        #smax     = sky + 2.5 * 1.5 * sigsky
        #smin     = sky - 2.5 * 1.0 * sigsky
        if diffimflag:
            smax     = sky + 2.5 * 1.5 * sigsky
            smin     = sky - 2.5 * 1.0 * sigsky
            scaling = imDisp.LimLinearScale(smin, smax)
        else:
            smax     = sky + 2.5 * 1.5 * sigsky
            smin     = sky - 2.5 * 1.0 * sigsky
            scaling = imDisp.LimLinearScale(smin, smax)

        # get header and scaling
        ##

        ##
        # get data limits and data
        x, y = list(map(float, sky2xy(fitsfile, ra, dec)))

        xmin = int(x - size/2)
        xmax = xmin + size

        ymin = int(y - size/2)
        ymax = ymin + size

        if xmin > xfits or ymin > yfits:
            # off the image totally!
            continue

        ffx0 = max(1, xmin)
        ffy0 = max(1, ymin)
        ffx1 = min(xfits, xmax)
        ffy1 = min(yfits, ymax)

        arx0 = ffx0 - xmin
        ary0 = ffy0 - ymin
        arx1 = size - (xmax - ffx1)
        ary1 = size - (ymax - ffy1)



        data = num.zeros((size, size))
        try:
            data[ary0:ary1, arx0:arx1] = fptr[0].data[ffy0:ffy1, ffx0:ffx1].copy()
        except:
            data = num.zeros((size, size))

        # get data limits and data
        ##


        if xflipflag:
            print('imagecutout: XFLIP')
            im    = Image.frombytes('P',(size, size), scaling(data).astype(num.int8).tobytes(),'raw','P').transpose(Image.FLIP_TOP_BOTTOM).transpose(Image.FLIP_LEFT_RIGHT)
        else:
            #im    = Image.fromstring('P',(size, size), scaling(data).astype(num.int8).tostring(),'raw','P').transpose(Image.FLIP_TOP_BOTTOM)
            im    = Image.frombytes('P',(size, size), scaling(data).astype(num.int8).tobytes(),'raw','P').transpose(Image.FLIP_TOP_BOTTOM)




        draw  = ImageDraw.Draw(im)
        draw.rectangle([0, 0, size, 11],   fill=128)
        #draw.text(     [5, 0], label, font=dfont, fill=255)
        color=255
        fcolor="cyan"
        if badimflag!=None:
            if badimflag[i]:
                color="red"
                fcolor=color
        #draw.text(     [5, 0], label, font=dfont, fill=color)

        ##draw.text(     [5, 0], flabel, font=dfont, fill=fcolor)
        ##draw.text(     [10, 0], label, font=dfont, fill=color)
        ##if diffimflag:
        ##    draw.text(     [5, 15], "DIFF", font=dfont, fill=fcolor)
        draw.text(     [5, 0], flabel, fill=fcolor)
        draw.text(     [10, 0], label,  fill=color)
        if diffimflag:
            draw.text(     [5, 15], "DIFF", fill=fcolor)

        if xhair:
            rad1 = 15
            rad2 = 25
            angle = math.sin(math.pi / 4.)

            #xc   = xmin + (xmax - xmin) / 2
            #yc   = ymin + (ymax - ymin) / 2
            xc   = size / 2
            yc   = size / 2

            for tx in [-1, 1]:
                for ty in [-1, 1]:
                    x1 = xc + tx * rad1 * angle
                    y1 = yc + ty * rad1 * angle
                    x2 = xc + tx * rad2 * angle
                    y2 = yc + ty * rad2 * angle
                    #draw.line( [x1, y1, x2, y2], fill=255)
                    #draw.line( [x1, y1, x2, y2], fill="red",width=4.0)
                    draw.line( [x1, y1, x2, y2], fill=fcolor,width=4)
        comp.paste(im, (X,Y))

        n += 1

        if n % nX == 0:
            X  = 0
            Y += size + buff1
        else:
            X += size + buff1

        del im, draw
        fptr.close()

    comp.save(outname)
    del font, comp, data

# Extend generic ConfigParser so that environment variables are substituted
class ConfigParser_env(SafeConfigParser):
    def __init__(self):
        SafeConfigParser.__init__(self)
        self.envvarpattern=re.compile('\$(\w+)')

    def getstring(self,section,paramname):
        s = self.subenvvarplaceholder(self.get(section,paramname))
        return(s)

    def gethex(self,section,paramname):
        val = int(eval(self.get(section,paramname)))
        return(val)

    def subenvvarplaceholder(self,s):
        envvarnames=self.envvarpattern.findall(s)
        if envvarnames:
            for name in envvarnames:
                envval=os.environ[name]
                subpattern='\$%s' % (name)
                s = re.sub(subpattern,envval,s)
        return(s)

    def setval_nosection(self,option,val,allflag=False,throwerror=True):
        sections = self.sections()
        foundflag=False
        errorflag=False
        for section in sections:
            if self.has_option(section,option):
                errorflag = self.set(section,option,val)
                if errorflag!=None and throwerror:
                    raise RuntimeError("ERROR %s %s %s" % (section,option,val))
                foundflag=True
                if not allflag:
                    break
        if (not foundflag) and throwerror:
            raise RuntimeError("ERROR: could not find section=%s, parameter=%s!" % (section,option))

        return(not foundflag)

    def setvals_nosection(self,option_val_list,allflag=False,throwerror=True):
        if option_val_list == None: return(0)
        errorflag = False
        for (option,val) in option_val_list:
            errorflag |= self.setval_nosection(option,val,allflag=allflag,throwerror=throwerror)
        return(errorflag)

    def setvals(self,section_option_val_list,throwerror=True):
        if section_option_val_list == None: return(0)
        errorflagall = False
        for (section,option,val) in section_option_val_list:
            if not self.has_option(section,option):
                errorflagall = True
                if throwerror:
                    raise RuntimeError("ERROR: section=%s, parameter=%s does not exist!" % (section,option))
                continue
            errorflag = self.set(section,option,val)
            if errorflag != None:
                errorflagall = True
                if throwerror:
                    raise RuntimeError("ERROR: could not set section=%s, parameter=%s to value %s!" % (section,option,val))
        return(errorflagall)

if __name__ == '__main__':
    ra  = '00:00:00.01'
    dec = '00:00:00.01'


