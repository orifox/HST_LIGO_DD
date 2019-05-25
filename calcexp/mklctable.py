#!/usr/bin/env python
import sys, os, math
# put the tools directory into the path
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
    tlc.loadfile('kilonova_phottable_40Mpc.txt')
    tlc.configcols(['f218w','f225w','f275w','f336w'],'f','%.1f')
    tlc.configcols(['rfphase','ofphase'],'f','%.2f')

    toutim = txttableclass()
    toutim.configcols(['phase'],'f','%.2f',visible=1)
    for filt in ['f218w','f225w','f275w','f336w']:
        toutim.configcols(['m(%s)' % filt],'f','%.1f',visible=1)
        toutim.configcols(['t(%s)' % filt],'f','%.0f',visible=1)

    #for specmode  in ['FUV,MAMA','NUV,MAMA','WFC3,G280']:
    #    toutim.configcols(['t(%s)' % specmode],'f','%.0f',visible=1)
        
    phases = [0.5,1.0,1.5,2.0,3.0,4.0,5.5,7.0]
    for phase in phases:
        toutim.newrow({'phase':phase})

    for filt in ['f218w','f225w','f275w','f336w']:
        tlc.initspline('rfphase',filt)
        for key in toutim.allrowkeys:
            m = tlc.spline(toutim.getentry(key,'phase'))
            toutim.setentry(key,'m(%s)' % filt,m)

    set_phase2t(toutim,1.0,zip(['f218w','f225w','f275w','f336w'],[700,150,30,5]))
    set_phase2t(toutim,1.5,zip(['f218w','f225w','f275w','f336w'],[3000,500,120,10]))
    set_phase2t(toutim,2.0,zip(['f218w','f225w','f275w','f336w'],[10000,2000,360,25]))
    set_phase2t(toutim,3.0,zip(['f225w','f275w','f336w'],[7000,1000,40]))
    set_phase2t(toutim,4.0,zip(['f225w','f275w','f336w'],[18000,2000,150]))
    set_phase2t(toutim,5.5,zip(['f275w','f336w'],[9000,300]))
    set_phase2t(toutim,7.0,zip(['f275w','f336w'],[15000,500]))


                       
    #set_phase2t(toutim,1.0,zip(['FUV,MAMA,','NUV,MAMA','WFC3,G280'],[t_1_orbit,1200,20*2]))
    #set_phase2t(toutim,1.5,zip(['NUV,MAMA','WFC3,G280'],[t_1_orbit,500]))
    #set_phase2t(toutim,2.0,zip(['WFC3,G280'],[2000]))
    #set_phase2t(toutim,3.0,zip(['WFC3,G280'],[t_1_orbit]))
    #set_phase2t(toutim,4.0,zip(['WFC3,G280'],[t_1_orbit]))
    
    toutim.printtxttable()

    toutspec = txttableclass()
    toutspec.configcols(['f275w'],'f','%.1f',visible=1)
    mags = [19.0,20.0,21.0,22.0,23.0,24.0]
    for m in mags:
        toutspec.newrow({'f275w':m})

    for instrument in ['STIS,CCD,G430L','STIS,CCD,G230LB','FUV,MAMA,G410L',
                       'FUV,MAMA,PRISM','NUV,MAMA,G230L','WFC3,G280']:
        toutspec.configcols([instrument],'f','%.0f',visible=1)
        
    
    set_m2t(toutspec,'STIS,CCD,G430L',zip([19.0,20.0,21.0,22.0],[160,420,1200,t_1_orbit]))
    set_m2t(toutspec,'STIS,CCD,G230LB',zip([19.0,20.0],[2600,10000]))
    set_m2t(toutspec,'FUV,MAMA,PRISM',zip([19.0,20.0,21.0,22.0,23.0,24.0],[15,40,120,320,1100,t_1_orbit]))
    set_m2t(toutspec,'WFC3,G280',zip([19.0,22.0,23.0,24.0],[20,500,2000,t_1_orbit]))
    set_m2t(toutspec,'FUV,MAMA,G410L',zip([19.0],[3*t_1_orbit]))
    set_m2t(toutspec,'NUV,MAMA,G230L',zip([19.0,20.0,21.0],[330,1200,t_1_orbit]))
    toutspec.printtxttable()

    for filt in ['f218w','f225w','f275w','f336w']:
        toutim.setcol2latexphantom('t(%s)' % filt)
    for instrument in ['STIS,CCD,G430L','STIS,CCD,G230LB','FUV,MAMA,G410L',
                 'FUV,MAMA,PRISM','NUV,MAMA,G230L','WFC3,G280']:
        toutspec.setcol2latexphantom(instrument)
    
    toutim.outputundefined='$\dots$'
    toutim.save2texfile('Im_Exptime.tex')
    toutspec.outputundefined='$\dots$'
    toutspec.save2texfile('Spec_Exptime.tex')

    
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
    tspecinfo.save2texfile('SpecInfo.tex')
