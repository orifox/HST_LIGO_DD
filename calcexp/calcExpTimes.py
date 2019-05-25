#!/usr/bin/env python
import sys, os, math, re, sys, copy
# put the tools directory into the path
#sys.path.append(os.path.join(os.environ['PIPE_PYTHONSCRIPTS'],'tools'))
#sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
from texttable import txttableclass
from tools import makepath4file
import argparse
import numpy as np

class calcexpclass(txttableclass):
    def __init__(self, *args, **kwargs):
        txttableclass.__init__(self, *args, **kwargs)

        self.verbose = False

        self.lc = txttableclass()
        self.lcfilters = []
        self.lc.magoffset = None

        self.phasecol = 'time'

        self.UVfilters = ['f218w','f225w','f275w','f336w']
        self.Ofilters = ['f438w','f606w','f814w']
        self.IRfilters = ['f110w','f160w']
        self.IRspec = ['WFC3_IR_G102','WFC3_IR_G141']
        self.NUVspec = ['WFC3_G280_2500','STIS_G230L_2500','STIS_PRISM_2500']
        self.FUVspec = ['STIS_PRISM_1500_SN3','STIS_G140L_1500_SN3']
        # which filter was used to calculate Spec exposure time
        self.filter4specmode = {'WFC3_IR_G102':'f110w','WFC3_IR_G141':'f110w','WFC3_G280_2500':'f275w','STIS_G230L_2500':'f275w','STIS_PRISM_2500':'f275w','STIS_PRISM_1500_SN3':'f218w','STIS_G140L_1500_SN3':'f218w'}

        self.UVoutcols = ['f218w','f225w','f275w','f336w']
        
        self.imaging_exptime = txttableclass()
        self.IRspec_exptime = txttableclass()
        self.UVspec_exptime = txttableclass()

        self.reftime = None
        self.reffilter = None
        self.refmag = None
        self.dt_KN = None
        self.dt_trigger = None

        self.lc_maxmag = 31.0
         
    def calc_magoffset(self,reftime,reffilter,refmag):
        self.reftime = reftime
        self.reffilter = reffilter
        self.refmag = refmag

        goodlckeys = self.lc.CUT_none_vals(reffilter)
        self.lc.initspline(self.phasecol,reffilter,keys=goodlckeys)
        m = self.lc.spline(reftime)
        self.lc.magoffset = refmag - m
        print('mag in filter %s at phase %.2f in lcfile %s: %.2f' % (reffilter,reftime,self.lc.filename,m))
        print('reference mag: %.2f, thus adding %.2f to all magnitudes of lc in lcfile' % (refmag,self.lc.magoffset))
        for col in self.lc.cols:
            if re.search('^f',col)==None:
                continue
            print(col)
            for key in self.lc.allrowkeys:
                if self.lc.getentry(key,col)==None:
                    continue
                #print(key,col,self.lc.getentry(key,col))
                self.lc.setentry(key,col,self.lc.getentry(key,col)+self.lc.magoffset)
        if self.verbose:
            self.lc.printtxttable()

        return(0)

    def loadlc(self,lcname):
        print('Loading %s' % lcname)
        self.lc.setpattern4undefined('^nan$')    
        if self.lc.loadfile(lcname):
            raise RuntimeError("Could not load %s" % lcname)
        self.lc.configcols(self.phasecol,'f','%.2f')
        self.lcfilters = []
        for col in self.lc.cols:
            if (re.search('^f',col)==None) and (re.search(col,'UBVRIugriz')==None) and (not(col in ['uvw2','uvm2','uvw1'])):
                continue
            self.lcfilters.append(col)
            self.lc.configcols(col,'f','%.2f')
            # some of the faint mags get unreliable, remove them
            keys = self.lc.CUT_inrange(col,self.lc_maxmag,None)
            for key in keys: self.lc.setentry(key,col,None)
        print('lcfilters: '+','.join(self.lcfilters))
        return(0)

    def load_exptime_table(self,tt,filename):
        tt.setpattern4undefined('^nan$')    
        print('Loading %s' % filename)
        if tt.loadfile(filename):
            raise RuntimeError("Could not load %s" % filename)
        if 'mag' in tt.cols: tt.configcols('mag','f','%.1f')
        if 'm_f275w' in tt.cols: tt.configcols('m_f275w','f','%.1f')
        if 'm_f218w' in tt.cols: tt.configcols('m_f218w','f','%.1f')
        for col in tt.cols:
            if col in ['mag','m_f275w','m_f218w','_id','_mask']:continue
            print col
            tt.configcols(col,'f','%.0f')

    def load_imaging_exptime_table(self,filename):
        self.load_exptime_table(self.imaging_exptime,filename)
        if self.verbose:
            self.imaging_exptime.printtxttable()

    def load_IRspec_exptime_table(self):
        #Annalisa put the IR spec exposure times into the imaging exposure table. Let's just link them then...
        self.IRspec_exptime = self.imaging_exptime
        
    def load_UVspec_exptime_table(self,filename):
        #Annalisa put the IR spec exposure times into the imaging exposure table. Let's just link them then...
        self.load_exptime_table(self.UVspec_exptime,filename)
        if self.verbose:
            self.UVspec_exptime.printtxttable()
        
    def init_outputtable(self,phasemin,phasemax,phasestep,dt_KN,dt_trigger):
        if dt_KN==None:
            self.dt_KN = self.reftime
        else:
            self.dt_KN = dt_KN

        self.dt_trigger = dt_trigger
        
        self.configcols(['phaseGW','phaseKN','phaseTrig'],'f','%.2f',visible=1)
        for phaseGW in np.arange(phasemin,phasemax,phasestep):
            phaseKN =  phaseGW - self.dt_KN
            phasetrigger = phaseKN - self.dt_trigger
            self.newrow({'phaseGW':phaseGW,'phaseKN':phaseKN,'phaseTrig':phasetrigger})
        #if self.verbose:
        #    self.printtxttable()

    def calcexptime_UVspec(self,specmodes):
        for specmode in specmodes:
            reffilter = self.filter4specmode[specmode]
            magcol = 'm_%s' % reffilter
            #print('VVVV',specmodes,magcol,self.UVspec_exptime.cols,specmode)
            self.configcols(['t_%s' % specmode],'f','%.0f',visible=1)

            # init light curve spline for given filter
            # get the good keys (i.e. no none vals) for the given light curve
            goodlckeys = self.lc.CUT_none_vals(reffilter)
            self.lc.initspline(self.phasecol,reffilter,keys=goodlckeys)
            maxmag4lc = self.lc.maxentry(reffilter,keys=goodlckeys)+0.1
            maxphase4lc = self.lc.maxentry(self.phasecol,keys=goodlckeys)

            # get the good keys (i.e. no none vals) for the given specmode in the exposure time table
            goodkeys = self.UVspec_exptime.CUT_none_vals(specmode)
            # init the spline in the exposure time table for the given specmode
            self.UVspec_exptime.initspline(magcol,specmode,keys=goodkeys,interpolationtype='logspline')

            # get the maximum magnitude for which we have an exposure time value for given specmode
            maxmag4exptime = np.amin((maxmag4lc,self.UVspec_exptime.maxentry(magcol,keys=goodkeys)+0.1))
            if self.verbose:
                print('Max mag for %s:' % specmode,maxmag4exptime)
            
            for key in self.allrowkeys:
                m = self.lc.spline(self.getentry(key,'phaseGW'))
                if m<=maxmag4exptime:
                    exptime = self.UVspec_exptime.spline(m)
                else:
                    exptime = None
                # never go beyond what we have in the light curve
                if self.getentry(key,'phaseGW')>maxphase4lc:
                    exptime = None
                #print(filt,self.getentry(key,'phaseGW'),m,exptime)
                self.setentry(key,'t_%s' % specmode,exptime)

    def calcexptime_IRspec(self):
        for specmode in self.IRspec:
            reffilter = self.filter4specmode[specmode]
            self.configcols(['t_%s' % specmode],'f','%.0f',visible=1)
            # init light curve spline for given filter
            goodlckeys = self.lc.CUT_none_vals(reffilter)
            self.lc.initspline(self.phasecol,reffilter,keys=goodlckeys)
            maxmag4lc = self.lc.maxentry(reffilter,keys=goodlckeys)+0.1
            maxphase4lc = self.lc.maxentry(self.phasecol,keys=goodlckeys)

            # get the good keys (i.e. no none vals) for the given specmode in the exposure time table
            goodkeys = self.IRspec_exptime.CUT_none_vals('t_%s' % specmode)

            # init the spline in the exposure time table for the given specmode
            self.IRspec_exptime.initspline('mag','t_%s' % specmode,keys=goodkeys,interpolationtype='logspline')

            # get the maximum magnitude for which we have an exposure time value for given specmode
            #maxmag4exptime = self.IRspec_exptime.maxentry('mag',keys=goodkeys)+0.1
            maxmag4exptime = np.amin((maxmag4lc,self.IRspec_exptime.maxentry('mag',keys=goodkeys)+0.1))
            if self.verbose:
                print('Max mag for t_%s:' % specmode,maxmag4exptime)
            
            #print 'VVV',specmode
            #if specmode == 'WFC3_IR_G102':
            #    print np.arange(18.0,23.0,0.1)
            #    for m in np.arange(18.0,23.0,0.1):
            #        print m,self.IRspec_exptime.spline(m)
            #    sys.exit(0)
                
            #print 'VVV',specmode,reffilter
            #if specmode == 'WFC3_IR_G102':
            #    print np.arange(0.0,10.0,0.11)
            #    for t in  np.arange(0.0,10.0,0.11):
            #        print t,self.lc.spline(t)
            #    sys.exit(0)
                
            for key in self.allrowkeys:
                m = self.lc.spline(self.getentry(key,'phaseGW'))
                if m<=maxmag4exptime:
                    exptime = self.IRspec_exptime.spline(m)
                else:
                    exptime = None
                # never go beyond what we have in the light curve
                if self.getentry(key,'phaseGW')>maxphase4lc:
                    exptime = None
                #print(filt,self.getentry(key,'phaseGW'),m,exptime)
                self.setentry(key,'t_%s' % specmode,exptime)
        
    def calcexptime_imaging(self):
        filters = copy.deepcopy(self.UVfilters)
        filters.extend(self.Ofilters)
        filters.extend(self.IRfilters)
        for filt in filters:
            # init new columns
            self.configcols(['m_%s' % filt],'f','%.2f',visible=1)
            self.configcols(['t_%s' % filt],'f','%.0f',visible=1)

            # init light curve spline for given filter
            goodlckeys = self.lc.CUT_none_vals(filt)
            self.lc.initspline(self.phasecol,filt,keys=goodlckeys)
            maxmag4lc = self.lc.maxentry(filt,keys=goodlckeys)+0.1
            maxphase4lc = self.lc.maxentry(self.phasecol,keys=goodlckeys)

            # get the good keys (i.e. no none vals) for the given filter in the exposure time table
            goodkeys = self.imaging_exptime.CUT_none_vals('t_%s' % filt)

            # init the spline in the exposure time table for the given filter
            self.imaging_exptime.initspline('mag','t_%s' % filt,keys=goodkeys,interpolationtype='logspline')

            # get the maximum magnitude for which we have an exposure time value for given filter
            #maxmag4exptime = self.imaging_exptime.maxentry('mag',keys=goodkeys)+0.1
            maxmag4exptime = np.amin((maxmag4lc,self.imaging_exptime.maxentry('mag',keys=goodkeys)+0.1))
            if self.verbose:
                print('Max mag for t_%s:' % filt,maxmag4exptime)

            for key in self.allrowkeys:
                #print(self.getentry(key,'phaseGW'),filt)
                m = self.lc.spline(self.getentry(key,'phaseGW'))
                if m<=maxmag4exptime:
                    exptime = self.imaging_exptime.spline(m)
                else:
                    exptime = None

                #print(self.getentry(key,'phaseGW'),maxphase4lc,m,maxmag4lc,maxmag4exptime)
                # never go beyond what we have in the light curve
                if self.getentry(key,'phaseGW')>maxphase4lc:
                    m = None
                    exptime = None
                #print(filt,self.getentry(key,'phaseGW'),m,exptime)
                self.setentry(key,'m_%s' % filt,m)
                self.setentry(key,'t_%s' % filt,exptime)


    def getoutputcols_imaging(self, filts, showmags=False):
        cols=[]
        for filt in filts:
            if showmags:
                cols.append('m_%s' % filt)
            cols.append('t_%s' % filt)
        return(cols)
        
    def getoutputcols(self, showmags=False, saveNUV=True, saveFUV=False, saveIR=False, saveO=False, skip_savespec=False):
        cols = ['phaseGW','phaseKN','phaseTrig']
        if saveNUV:
            cols.extend(self.getoutputcols_imaging(self.UVfilters,showmags=showmags))
            if not skip_savespec:
                for specmode in self.NUVspec:
                    cols.append('t_%s' % specmode)
        if saveFUV and (not skip_savespec):
            for specmode in self.FUVspec:
                cols.append('t_%s' % specmode)
        if saveO:
            cols.extend(self.getoutputcols_imaging(self.Ofilters,showmags=showmags))
        if saveIR:
            cols.extend(self.getoutputcols_imaging(self.IRfilters,showmags=showmags))
            if not skip_savespec:
                for specmode in self.IRspec:
                    cols.append('t_%s' % specmode)
                        
        return(cols)
                    
        
    def savetable(self,cols=None):
        if cols is None:
            print('WARNING: No columns to save!! returning...')
            return(1)
        makepath4file(self.outbasename)
        print('Saving %s' % '%s.txt' % self.outbasename)
        self.save2file('%s.txt' % self.outbasename,cols=cols)

        return(0)
    
    def savelatextable(self, tt, outfilename=None, cols=None):
        if cols is None:
            print('WARNING: No columns to save!! returning...')
            return(1)

        if outfilename==None:
            outfilename = '%s.tex' % self.outbasename
            
        makepath4file(outfilename)
        tt.outputundefined='$\dots$'   
        for col in cols:
            if re.search('^t_',col) or (col in self.FUVspec) or (col in self.NUVspec) or (col in self.IRspec) :
                tt.setcol2latexphantom(col,latexphantomflag=True)

        print('Saving %s' % outfilename)
        tt.save2plaintextable(outfilename,cols=cols)

        tt.outputundefined='-'
        for col in cols:
            tt.setcol2latexphantom(col,latexphantomflag=False)
                
        return(0)
            
    def setoutbasename(self,outrootdir=None,outsubdir=None,outbasename=None,addsuffix=None,skip_refinfo=False):
        if outrootdir == None:
            outrootdir = '.'
        s = os.path.abspath(outrootdir)
        if outsubdir is not None:
            s += '/%s' % outsubdir
        if outbasename is None:
            (outbasename,suffix) = os.path.splitext(os.path.basename(self.lc.filename))
        s += '/%s' % outbasename

        if not skip_refinfo:
            s+='_t%.2f_%s_%.2f' % (self.reftime,self.reffilter,self.refmag)

        if addsuffix is not None:
            s += '.%s' % addsuffix

        self.outbasename = s
        print('outbasename: %s' % self.outbasename)
        return(0)
        
if __name__ == '__main__':
    calcExp = calcexpclass()

    parser = argparse.ArgumentParser()

    parser.add_argument("reftime", type=float, help=("time difference between KN and to GW discovery in days"))
    parser.add_argument("reffilter",  help=("reference filter"))
    parser.add_argument("refmag", type=float, help=("reference mag in reference filter"))

    parser.add_argument('-v', '--verbose', help="verbose level (default=%(default)s)",
                        action="store_true", default=False)
    parser.add_argument('-m', '--showmags', help="show the mags in the output tex table (default=%(default)s)",
                        action="store_true", default=False)
    parser.add_argument( '--savetable', help="save the table. columns depend on --saveNUV, --saveIR etc",
                        action="store_true", default=False)
    parser.add_argument( '--savelatex', help="save the latex table",
                        action="store_true", default=False)
    parser.add_argument( '--saveNUV', help="save NUV info into output",
                        action="store_true", default=False)
    parser.add_argument( '--saveFUV', help="save FUV info into output",
                        action="store_true", default=False)
    parser.add_argument( '--saveIR', help="save IR info into output",
                        action="store_true", default=False)
    parser.add_argument( '--saveO', help="save optical info into output",
                        action="store_true", default=False)
    parser.add_argument( '--skip_savespec', help="Skip the spec info in the output",
                         action="store_true", default=False)

    #parser.add_argument('--lc', default='./GW170817_40Mpc.txt',
    parser.add_argument('--lc', default='./kilonova_phottable_gw170817_40Mpc.txt',
                        help='lightcurve filename (default=%(default)s)')
    parser.add_argument('--imaging_IRspec_exptime', default='./exp_times_imaging_IRgrism.txt',
                        help='exposure time table for imaging and NIR Grism spectroscopy (default=%(default)s)')
    parser.add_argument('--UVspec_exptime', default='./exp_times_UVspec.txt',
                        help='exposure time table for UV spectroscopy (default=%(default)s)')

    parser.add_argument("--addoutsuffix", default=None, help=("time difference between KN and GW discovery in days. If not specified, then it is assumed that dt_KN=reftime"))

    parser.add_argument("--dt_KN", type=float, default=None, help=("time difference between KN and GW discovery in days. If not specified, then it is assumed that dt_KN=reftime"))
    parser.add_argument("--dt_trigger", type=float, default=0.1, help=("time difference between HST trigger and KN discovery in days"))
    parser.add_argument("--phasemin", type=float, default=0.0, help=("minimum phase since GW discovery in days (default=%(default)s)"))
    parser.add_argument("--phasemax", type=float, default=14.1, help=("maximum phase since GW discovery in days (default=%(default)s)"))
    parser.add_argument("--phasestep", type=float, default=0.5, help=("phase stepsize in days (default=%(default)s)"))

    parser.add_argument('--outrootdir', default=None,
                        help=('output root directory.'
                              '(default=%(default)s)'))
    parser.add_argument('--outbasename', default=None,
                        help=('output basename. If None, then the basename of input lc is taken'
                              '(default=%(default)s)'))
    parser.add_argument('--outsubdir', default=None,
                        help=('subdir added to the output root directory (and filename) '
                              '(default=%(default)s)'))
    parser.add_argument('--addsuffix', default=None,
                        help='suffix added to the output basename (default=%(default)s)')
    parser.add_argument('--skip_refinfo', help="Skip putting the ref info (reftime, refmag, reffilter) into the filename (default=%(default)s)",
                        action="store_true", default=False)


    
    args = parser.parse_args()
    calcExp.verbose=args.verbose

    # load lightcurve table and calculate delta mag to move light curves to reference mag in reference filter
    calcExp.loadlc(args.lc)
    calcExp.calc_magoffset(args.reftime,args.reffilter,args.refmag)

    calcExp.setoutbasename(outrootdir=args.outrootdir,outsubdir=args.outsubdir,
                           outbasename=args.outbasename,addsuffix=args.addsuffix,
                           skip_refinfo=args.skip_refinfo)
    # load imaging exposure time tables
    calcExp.load_imaging_exptime_table(args.imaging_IRspec_exptime)
    # load IRspec exposure time tables
    calcExp.load_IRspec_exptime_table()
    # load UVspec exposure time tables
    calcExp.load_UVspec_exptime_table(args.UVspec_exptime)

    if 1==1:
        outputdir=os.path.dirname(calcExp.outbasename)
        calcExp.savelatextable(calcExp.imaging_exptime,outfilename='%s/UVim_exptime.tex' % outputdir,cols=['mag','t_f218w','t_f225w','t_f275w','t_f336w'])
        calcExp.savelatextable(calcExp.IRspec_exptime,outfilename='%s/IR_exptime.tex' % outputdir,cols=['mag','t_f110w','t_f160w','t_WFC3_IR_G102','t_WFC3_IR_G141'])
        calcExp.savelatextable(calcExp.UVspec_exptime,outfilename='%s/UV_exptime.tex' % outputdir,cols=['m_f275w','WFC3_G280_2500','STIS_G230L_2500','STIS_PRISM_2500','m_f218w','STIS_PRISM_1500_SN3','STIS_G140L_1500_SN3'])
        calcExp.savelatextable(calcExp.UVspec_exptime,outfilename='%s/NUV_exptime.tex' % outputdir,cols=['m_f275w','WFC3_G280_2500','STIS_G230L_2500','STIS_PRISM_2500'])
        calcExp.savelatextable(calcExp.UVspec_exptime,outfilename='%s/FUV_exptime.tex' % outputdir,cols=['m_f218w','STIS_PRISM_1500_SN3','STIS_G140L_1500_SN3'])
        
    #initialize the output table: calculate the phases
    calcExp.init_outputtable(args.phasemin,args.phasemax,args.phasestep,args.dt_KN,args.dt_trigger)

    # fill up output table with  imaging exposure times and mags
    calcExp.calcexptime_imaging()

    # fill up output table with  IR spec exposure times 
    calcExp.calcexptime_IRspec()

    # fill up output table with  UV spec exposure times 
    calcExp.calcexptime_UVspec(calcExp.NUVspec)
    calcExp.calcexptime_UVspec(calcExp.FUVspec)

    outcols = calcExp.getoutputcols(showmags=args.showmags, saveNUV=args.saveNUV, saveFUV=args.saveFUV, saveIR=args.saveIR, saveO=args.saveO, skip_savespec=args.skip_savespec)
    
    if calcExp.verbose:
        calcExp.printtxttable(cols=outcols)

    # Save the table
    if args.savetable:
        calcExp.savetable(cols=outcols)

    # Save the latextable
    if args.savelatex:
        calcExp.savelatextable(calcExp,cols=outcols)
