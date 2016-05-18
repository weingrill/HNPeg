#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 10, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
from datasource import DataSource

def getfile(sourcefile, targetdirectory):
    from subprocess import call
    call(['scp', 'sro@papaya.aip.de:'+sourcefile,targetdirectory])
    pass

def convertfile(filename):
    from subprocess import call
    call(['./ses-writetxt.py', '--textfile', filename])

if __name__ == '__main__':
    query = "SELECT filename FROM obs WHERE object LIKE 'HN Peg' ORDER BY dateobs;"
    wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
    result = wifsip.query(query)
    for r in result: 
        filename = r[0]
        print filename,
        #science20151106B-0038_botzfxsEcd.fits
        path = filename.lstrip('science')[:8]
        sourcefile = '/stella/home/stella/spectra/'+path+'/'+filename+'_botzfxsEcd.fits'
        targetdirectory = '/work2/jwe/Projects/HNPeg/data'
        print sourcefile, targetdirectory
        #getfile(sourcefile,targetdirectory)
        convertfile('/work2/jwe/Projects/HNPeg/data/'+filename+'_botzfxsEcd.fits')
        #/stella/home/stella/spectra/20151106
        

    print len(result)