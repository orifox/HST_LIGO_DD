#!/usr/bin/env python
import sys, os, math
# put the tools directory into the path
sys.path.append(os.path.join(os.environ['PIPE_PYTHONSCRIPTS'],'tools'))
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
from texttable import txttableclass

def set_phase2t(t,phase,filtlist):
    keys = t.CUT_inrange('phase',phase,phase)
    if len(keys)!=1:
        raise RuntimeError("Could not find phase %s!" % phase)
    for (filt,time) in filtlist:
        t.setentry(keys[0],'t(%s)' % filt,time)
    return(0)
        
def set_m2t(t,instrument,m2t):
    for (m,time) in m2t:
        keys = t.CUT_inrange('f275w',m,m)
        if len(keys)!=1:
            raise RuntimeError("Could not find m %f!" % m)
        t.setentry(keys[0],'%s' % instrument,time)
    return(0)
def set_specinfo(t,instrument,dict):
    for (row,info) in dict:
        keys = t.CUT_inrange('X',row,row)
        if len(keys)!=1:
            raise RuntimeError("Could not find row %s!" % row)
        t.setentry(keys[0],'%s' % instrument,info)
    return(0)
        
if __name__ == '__main__':

    t_1_orbit = 3000.0
    
    tlc = txttableclass()
    tlc.loadfile('kilonova_phottable_40Mpc_v2.txt')
    tlc.configcols(['f218w','f225w','f275w','f336w'],'f','%.1f')
    tlc.configcols(['rfphase','ofphase'],'f','%.2f')

    toutim = txttableclass()
    toutim.setpattern4undefined('^nan$')    
    toutim.loadfile('exptimes_UV_imaging.txt')
    toutim.configcols(['t_f218w','t_f225w','t_f275w','t_f336w','t_t25cn182','t_t25cn270'],'d','%d',latexphantomflag=True)
    toutim.configcols(['mag'],'f','%.0f',visible=1)
    #toutim.printtxttable()


    #for filt in ['t_f218w','t_f225w','t_f275w','t_f336w','t_t25cn182','t_t25cn270']:
        #toutim.colinfo[filt]['latexphantomflag']=True
    #    toutim.setcol2latexphantom(filt,latexphantomflag=True)
    toutim.outputundefined='$\dots$'

    print('Saving UVimaging.tex')
    toutim.save2plaintextable('UVimaging.tex')

    # undo output changes for latex...
    toutim.outputundefined='-'
    for filt in ['t_f218w','t_f225w','t_f275w','t_f336w','t_t25cn182','t_t25cn270']:
        toutim.setcol2latexphantom(filt,latexphantomflag=False)

    toutim.printtxttable()

    tspecinfo = txttableclass()
    tspecinfo.configcols(['X','STIS,CCD,G430L','STIS,CCD,G230LB','FUV,MAMA,G410L',
                          'FUV,MAMA,PRISM','NUV,MAMA,G230L','WFC3,G280'],'s','%s',visible=1)

    specinforows = ['Instrument','Detector','Grating','Type','Central $\lambda$','Range','Slit','Resolution','F275W','Comments']
    for row in specinforows:
        tspecinfo.newrow({'X':row})
        
    set_specinfo(tspecinfo,'STIS,CCD,G430L',zip(specinforows,
                                                ['STIS','CCD','G430L','NUV','FILLME','2900 - 5700','52x2"','530 - 1040','22.0']))
    set_specinfo(tspecinfo,'STIS,CCD,G230LB',zip(specinforows,
                                                ['STIS','CCD','G280','NUV','FILLME','1685 - 3065','52x2"','620 - 1100','20.0']))
    set_specinfo(tspecinfo,'FUV,MAMA,G410L',zip(specinforows,
                                                ['STIS','MAMA','G410L','FUV','FILLME','1150 - 1736','52x2"','960 - 1440','18.0','low glow']))
    set_specinfo(tspecinfo,'FUV,MAMA,PRISM',zip(specinforows,
                                                ['STIS','MAMA','PRISM','FUV','FILLME','1150 - 3000','slitless"','10 - 2500','24.0','low glow']))
    set_specinfo(tspecinfo,'NUV,MAMA,G230L',zip(specinforows,
                                                ['STIS','MAMA','G230L','NUV','FILLME','1600 - 3000','52x2"','500 - 1010','19.0']))
    set_specinfo(tspecinfo,'WFC3,G280',zip(specinforows,
                                                ['WFC3','UVIS','G280','NUV','FILLME','1900 - 4500','slitless','70 - 300','24.0']))
    
    tspecinfo.printtxttable()
    tspecinfo.save2plaintextable('SpecInfo.tex')
