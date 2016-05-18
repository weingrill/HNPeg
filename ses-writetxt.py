#!/usr/bin/env python

import pyfits
import numpy
import re
import scipy.special
import scipy.interpolate
import optparse

# speed of light in km/s
C = 299792.458

# SES order to aperture (starting at 1) offset: add 1 for c-style array iteration
offset = 64 #(order 1 is order 65)

def join_struct_arrays(arrays):
    newdtype = sum((a.dtype.descr for a in arrays), [])
    newrecarray = numpy.empty(len(arrays[0]), dtype = newdtype)
    for a in arrays:
        for name in a.dtype.names:
            newrecarray[name] = a[name]
    return newrecarray

def join_struct_arrays2(arrays):
    sizes = numpy.array([a.itemsize for a in arrays])
    offsets = numpy.r_[0, sizes.cumsum()]
    n = len(arrays[0])
    print sizes, offsets, n
    joint = numpy.empty((n, offsets[-1]), dtype=numpy.float64)
    for a, size, offset in zip(arrays, sizes, offsets):
        joint[:,offset:offset+size] = a.view(numpy.float64).reshape(n,size)
    dtype = sum((a.dtype.descr for a in arrays), [])
    return joint.ravel().view(dtype)


def calcwavelength(p,npts,z):
  #     w = sum from i=1 to nfunc {wt_i * (w0_i + W_i(p)) / (1 + z)}
  # we allocate the wavelength vector
  ww = numpy.zeros( (npts) )
  xn = numpy.arange(int(p[0][4]),int(p[0][5]+1), dtype=int)
  x = (2.*xn-(p[0][5]+p[0][4]))/(p[0][5]-p[0][4])
  # we create the needed poynomial functions
  Pn = []
  for i in range(p[0][3]):
    if p[0][2] == 1:
      # chebychev polynomials
      Pn.append( scipy.special.chebyt(i))
    else:
      print "unsopported dispersion function"
      return None
  for n in range(len(p)):
    spec_w0 = p[n][1]
    for i in range(p[n][3]):
      ww += p[n][0]* (spec_w0 + p[n][i+6]*Pn[i](x)) / (1.+z)
  return ww


# Function that returns an array of wavelength
# hdr is a pyfits FITS-header 
# od  is the corresponding data section (used for getting the dimensions)
# returns an array the size of od with wavelengths for each pixel
# NOTE: not all dispersion possibilities of IRAF are implemented
#       as of June 2007: echelle with linear dispersion and 4th order cheb. polynomials
def GetWavelength(hdr, od, extended=False):
  # multidimensional is probably Echelle
  if od.ndim > 1:
    (norder, npts) = od.shape
    # lets see if it is OHP format
    try:
      lloffset=float(hdr['LLOFFSET'])
      degxll = int(hdr['DEGXLL'])+1
      degoll = int(hdr['DEGOLL'])+1
      llpars = float(hdr['LLPARS'])
      coell = numpy.zeros( (degxll, degoll) )
      coell[0,0] = 0.25 * float(hdr['COELL1'])
      for i in range(1,degxll):
        coell[i,0] = 0.5 * float(hdr['COELL%d' % (1+i*(degoll))])
      for j in range(1,degoll):
        coell[0,j] = 0.5 * float(hdr['COELL%d' % (1+j)])
      for i in range(1,degxll):
        for j in range(1,degoll):
          coell[i,j] = float(hdr['COELL%d' % (1+j+i*(degoll))])
      m = numpy.zeros(od.shape)
      x = numpy.zeros(od.shape)
      for i in range(norder):
        m[i] = lloffset + llpars * float(i+1)
      for j in range(npts):
        x[:,j] = float(j+1)
      allapertures = numpy.int_(lloffset + llpars * (numpy.arange(norder)+1))
      dw = dispersionf_test(coell,m,x,(degxll, degoll))
    except:
      # No? We assume it is a iraf WAT format header
      wat = ""
      for i in range(1,1000):
        try:
          wat += hdr["WAT2_%03d" % (i)].ljust(68)
        except KeyError:
          break
      watarr = re.findall('("[^"]+")+', wat)
      dw = numpy.zeros([norder, npts], numpy.float64)
      allapertures = numpy.zeros( (norder), int)
      for ord in range(norder):
        #  first we check what happens if we assume a one digit exponential
        nn = re.findall('[-+]?[0-9]*\.?[0-9]+[eE][-+]?[0-9]', watarr[ord])
        longexp = False
        for i in range(len(nn)):
          # if there is at least one bogus exponential, we allow longer exponentials 
          if nn[i][-1] == "0":
            longexp = True
            break
        if not longexp:
          # this is iraf-bug-safe, but crashes if exponentials have more than 1 digit, i.e. "1.e-01" instead of "1.e-1"
          nn = re.findall('[-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]?', watarr[ord])
        else:
          nn = re.findall('[-+]?[0-9]*\.?[0-9]+[eE]?[-+]?[0-9]*', watarr[ord])
        #print len(nn), nn
        order = int(nn[0])        # order number, is same as aperture for old setup
        spec_beam = int(nn[1])    # beam, unused
        allapertures[ord] = spec_beam
        #spec_offset = spec_beam - ord - 1  # offset between aperture and order number (not beam!)
        type = int(nn[2])    # 0=linear, 1= log-linear,2=nonlinear
        spec_w1 = float(nn[3])
        spec_dw = float(nn[4])
        spec_nw = int(nn[5])
        spec_z = float(nn[6])
        if len(nn)>7:
          spec_aplow = float(nn[7])
          spec_aphigh = float(nn[8])
        else:
          spec_aplow = 0.
          spec_aphigh = 0.
        #w = (w1 + dw * (p - 1)) / (1 + z)
        #w = scipy.special.exp10 {(w1 + dw * (p - 1)) / (1 + z)}
        if type == 0:
	  for l in range(npts):
	    dw[ord,l] = (spec_w1 + float(l) * spec_dw) / (1.+spec_z)
        elif type == 1:
	  for l in range(npts):
	    dw[ord,l] = 10.**( (spec_w1 + float(l) * spec_dw) / (1.+spec_z))
        elif type == 2:
          # nonlinear
          nfunc = 0
          ncur = 9
          spec_funcs = []
          while ncur < len(nn):
            spec_funcs.append([])
            spec_funcs[nfunc].append(float(nn[ncur]))    # weight
            spec_funcs[nfunc].append(float(nn[ncur+1]))    # w0
            spec_funcs[nfunc].append(int(nn[ncur+2]))    # type
            polyorder = int(nn[ncur+3])
            #look for E-10 = E-1 + 0 style problems
            if len(nn) > ncur+6+polyorder:
              if ( nn[ncur+6+polyorder]== '0'):
                nn[ncur+6+polyorder-1]+='0'
                del(nn[ncur+6+polyorder])
            spec_funcs[nfunc].append(polyorder)    # order
            spec_funcs[nfunc].append(float(nn[ncur+4]))    # xmin
            spec_funcs[nfunc].append(float(nn[ncur+5]))    # xmax
            for i in range(polyorder):
              spec_funcs[nfunc].append(float(nn[ncur+6+i]))
            ncur += 6 + polyorder
            nfunc += 1
          dw[ord,:] = calcwavelength(spec_funcs, npts, spec_z)
        else:
          print "Error, unsupported data format: %d" % type
  else:
       # 1d spectra are always linear?
       #spec_offset = 0
       npts = od.size
       lam0_t = float(hdr["CRVAL1"])
       dlam_t = float(hdr["CDELT1"])
       try:
         ltv = float(hdr["LTV1"])
       except:
         ltv=0
       dw = numpy.empty((npts))
       for l in range(npts):
           dw[l] = lam0_t + float(l-ltv) * dlam_t
  if extended == True:
    return dw, allapertures
  else:
    return dw


from optparse import OptionParser
usage = "usage: %prog [options] ses-fits-files"
parser = OptionParser(usage=usage, version="%prog 1.0")   
parser.add_option("-v", "--vcorr",
                  action="store", dest="vcorr", type="int", 
                  default=3, help="Apply radial velocity shift: 0=none, 1=vobs, 2=vorbit, 3=vhelio_corr. default:3")
parser.add_option("--fitstable", action="store_true", dest="fitstable", default=False, help="Create FITS table instead of text file")
parser.add_option("--fitsfile", action="store_true", dest="fitsfile", default=False, help="Create log-linear sampled FITS file instead of text file")
parser.add_option("--textfile", action="store_true", dest="textfile", default=False, help="Create log-linear sampled text file instead of standard text file")
parser.add_option("--oversample", action="store", dest="oversample", default=2, type="int", help="Oversampling for equally sampled spectrum")

(opts, args) = parser.parse_args()
if len(args) < 1:
  parser.error("Incorrect number of arguments")


for n in range(len(args)):
  fname = args[n]
  try:
    #ffile = pyfits.open(fname, mode='update')
    ffile = pyfits.open(fname)
  except:
    print "Error opening %s" % fname
  else:
    hdr = ffile[0].header;
    if ffile[0].data.ndim == 1:
      nblock,norder,npts = 1, 1, ffile[0].data.shape[0] 
      od = numpy.empty( (1,npts))
      od[0]  = ffile[0].data
      od_cont  = od.copy()
      Err  = numpy.zeros(od.shape) + 0.01
      lam0 = float(hdr["CRVAL1"])
      try:
        dlam = float(hdr["CDELT1"])
      except:
        dlam = float(hdr["CD1_1"])
      try:
        ltv = float(hdr["LTV1"])
      except:
        ltv=0
      try:
        crpix = float(hdr["CRPIX1"])
      except:
        crpix=1.
      wd = numpy.empty( (1,npts))
      wd[0] = numpy.arange(lam0-(crpix-1.)*dlam, lam0-(crpix-1.)*dlam+(npts-0.5)*dlam, dlam)
    else:
      (nblock,norder,npts) = ffile[0].data.shape
      if nblock == 2:
        nod = 0
        nodc = 0
        nerr = 1
        od  = ffile[0].data[0]
        od_cont  = ffile[0].data[0]
        Err  = od/ffile[0].data[1]
      else:
        od  = ffile[0].data[0]
        od_cont  = ffile[0].data[1]
        Err  = ffile[0].data[2]
      wd  = GetWavelength(hdr,od) 
    ffile.close()
    vshift = 0.
    ## we get vhelio from the first spectrum to have a reference for vcorr=2 
    #if n == 0:
    #   try:
    #      vhelio = float(hdr['VHELIO'])
    #   except:
    #      vhelio = 0.0
    vhelio = 0.0

    if opts.vcorr == 1:
       try:
          vshift = float(hdr['VOBS'])
       except:
          vshift = 0.0
    elif opts.vcorr == 2:
       vshift = float(hdr['VORBIT'])
    elif opts.vcorr == 3:
       vshift = vhelio-float(hdr['VCORRECT'])
    if vshift != 0.:
       wd = wd / (1. + vshift/C)
    #print "correct for v = %s" % (vshift)
    (norder, npts) = od.shape
    sn = od / Err
    # we make sure that we do not divide by zero
    # cont_cond is later also used for the interpolating spline fit 
    cont_cond = od_cont > numpy.finfo(float).eps
    cont_ratio = numpy.empty(od_cont.shape)
    cont_ratio.fill(numpy.nan)
    cont_ratio[cont_cond] = od[cont_cond] / od_cont[cont_cond]
    err_cont = numpy.empty(od_cont.shape)
    err_cont.fill(numpy.nan)
    err_cont[cont_cond] = Err[cont_cond] / cont_ratio[cont_cond]
    try:
      trimarr_file = pyfits.open('trimarr%dx%d.fits' % (norder,npts))
    except:
      trimarr_file = pyfits.open('trimarr.fits')
    trimarr=trimarr_file[0].data
    trimhdr=trimarr_file[0].header
    t = (trimarr > 0)
    t2 = (trimarr > 0) & cont_cond
    trimarr_file.close()
    norder = trimhdr['NORDER']
    nslices = trimhdr['NSLICES']
    
    if opts.fitstable or opts.fitsfile or opts.textfile:
      rebins = numpy.zeros((norder/nslices), dtype=int)
      rebine = numpy.zeros((norder/nslices), dtype=int)
      for i in range(norder/nslices):
        rebins[i] = trimhdr['REBINS%d' % i]
        rebine[i] = trimhdr['REBINE%d' % i]
      '''
      # we create an array with the trimmed data to sort for wavelengths
      wave = numpy.array((), dtype=float)
      int_norm = numpy.array((), dtype=float)
      err_norm = numpy.array((), dtype=float)
      for i in range(norder):
        wave = numpy.append(wave, wd[i][t[i]])
        int_norm = numpy.append(int_norm, od_cont[i][t[i]])
        err_norm = numpy.append(err_norm, err_cont[i][t[i]])
      specarr = numpy.core.records.fromarrays(numpy.vstack((wave, int_norm, err_norm)), 
                                             names='x, y, e',
                                             formats = 'f8, f8, f8')
      specarr.sort(order='x')
      '''
      if opts.fitsfile or opts.textfile:
        # the fitsfile is linearly log-lambda spaced
        steps =  numpy.log10(wd[:,1:]/numpy.roll(wd,1,axis=1)[:,1:])
        step = numpy.median(steps)/float(opts.oversample)
        nptsNew = int(numpy.log10(wd[t].max()/wd[t].min())/step)
        newx = wd[t].min() * 10.**(step*numpy.arange(nptsNew, dtype=float))
        #print nptsNew, newx[226355-1], newx.max(),  wd[t].max()
        #nknots = int(numpy.log10(specarr['x'].max()/specarr['x'].min())/(step*3.))
        #knots = specarr['x'].min() * 10.**(3.*step*numpy.arange(nknots, dtype=float))
        #s = scipy.interpolate.LSQUnivariateSpline(specarr['x'],specarr['y'], t=knots[1:-2],w=1./specarr['e'])
        #nptsNew = int(numpy.log10(specarr['x'].max()/specarr['x'].min())/step)
        #newx = specarr['x'].min() * 10.**(step*numpy.arange(nptsNew, dtype=float))
      else:
        # the fitstable is made from wavelengths at original positions
        newx = numpy.empty((0))
        for i in range(norder/nslices-1, -1, -1):
          ap = i*nslices
          newx = numpy.append(newx, wd[ap][rebins[i]:rebine[i]])
      #print len(newx)
      # construct spline interpolators
      splinterpol = []
      splinterpol_err = []
      # construct structured array with the edges of the orders
      wave1 = numpy.zeros((norder*2))
      typ1 = numpy.zeros((norder*2), dtype=int)
      ordernum = numpy.zeros((norder*2), dtype=int)
      spos = numpy.zeros((norder*2), dtype=int)
      for i in range(norder):
        wave1[i] = wd[i][t[i]][0]
        typ1[i] = 0
        ordernum[i] = i
        wave1[norder+i] = wd[i][t[i]][-1]
        typ1[norder+i] = 1
        ordernum[norder+i] = i
        spos[i] = numpy.searchsorted(newx, wave1[i])
        spos[norder+i] = numpy.searchsorted(newx, wave1[norder+i])
        #spos[norder+i] = numpy.searchsorted(newx, wave1[norder+i])-1
        #print wave1[i], newx[spos[i]], wave1[norder+i], newx[spos[norder+i]]
        splinterpol.append( scipy.interpolate.InterpolatedUnivariateSpline(wd[i][t2[i]],od_cont[i][t2[i]]))
        splinterpol_err.append( scipy.interpolate.InterpolatedUnivariateSpline(wd[i][t2[i]],err_cont[i][t2[i]]))
        #splinterpol_err.append( scipy.interpolate.InterpolatedUnivariateSpline(wd[i][t[i]],Err[i][t[i]]))

      orderarr = numpy.core.records.fromarrays(numpy.vstack((wave1, typ1, ordernum, spos)), 
                                             names='x, t, o, p',
                                             formats = 'f8, i4, i4, i4')
      orderarr.sort(order='x')
      # step through the pieces between order boundaries and average the spline-interpolated orders
      newy = numpy.empty(newx.shape)
      newy.fill(numpy.nan)
      newerr = numpy.empty(newx.shape)
      newerr.fill(numpy.nan)
      curlist = []
      for i in range(len(orderarr)-1):
        if orderarr[i]['t'] == 0: 
          curlist.append(orderarr[i]['o'])
        else:
          curlist.remove(orderarr[i]['o'])
        totalweight = numpy.zeros(newx[orderarr[i]['p']:orderarr[i+1]['p']].shape)
        newy[orderarr[i]['p']:orderarr[i+1]['p']] = 0.
        newerr[orderarr[i]['p']:orderarr[i+1]['p']] = 0.
        for j in range(len(curlist)):
          newy[orderarr[i]['p']:orderarr[i+1]['p']] += splinterpol[curlist[j]](newx[orderarr[i]['p']:orderarr[i+1]['p']]) / splinterpol_err[curlist[j]](newx[orderarr[i]['p']:orderarr[i+1]['p']])
          newerr[orderarr[i]['p']:orderarr[i+1]['p']] += splinterpol_err[curlist[j]](newx[orderarr[i]['p']:orderarr[i+1]['p']])**2
          totalweight += 1./splinterpol_err[curlist[j]](newx[orderarr[i]['p']:orderarr[i+1]['p']])
        if len(curlist) > 0:
          newy[orderarr[i]['p']:orderarr[i+1]['p']] /= totalweight
          #newerr[orderarr[i]['p']:orderarr[i+1]['p']] = numpy.sqrt(newerr[orderarr[i]['p']:orderarr[i+1]['p']]) / float(len(curlist))
          #print newx[orderarr[i]['p']:orderarr[i+1]['p']], newy[orderarr[i]['p']:orderarr[i+1]['p']]
          #print orderarr[i]['p'],orderarr[i+1]['p']
      #print newy[-2:]
      if opts.textfile:
        outfile = open(fname.replace(".fits", ".txt"), "w+")
        try:
          outfile.write("# OBJNAME = %s\n" % (hdr['OBJNAME']))
        except:
          outfile.write("# OBJNAME = %s\n" % (hdr['OBJECT']))
        try:
          outfile.write("# HJD = %s\n" % (hdr['HJD']))
        except:
          pass
        try:
          outfile.write("# VOBS = %s\n" % (hdr['VOBS']))
        except:
          pass
        try:
          outfile.write("# VORBIT = %s\n" % (hdr['VORBIT']))
        except:
          pass
        try:
          outfile.write("# VCORRECT = %s\n" % (hdr['VCORRECT']))
        except:
          outfile.write("# VCORRECT = %s\n" % (float(hdr['VHELIO'])-float(hdr['VOBS'])))
        for n in range(len(newx)):
            outfile.write("%.4f %.3f %.3f\n" % (newx[n], newy[n], newerr[n]))
        outfile.close()
      elif opts.fitsfile:
        #cond = numpy.isnan(newy)
        #print len(numpy.compress(cond, newy))
        #print newx,newy, newerr
        hdu_out = pyfits.PrimaryHDU(numpy.vstack((newy, newerr)))
        hdr_out = hdu_out.header
        hdr_out.extend( hdr.copy(strip=True))
        for i in range(1,1000):
          if hdr_out.__contains__("WAT2_%03d" % (i)):
            del hdr_out["WAT2_%03d" % (i)]
          else:
            break
        hdr_out['DC-FLAG'] = 1
        hdr_out['APNUM1'] = '1 1' 
        hdr_out['APNUM2'] = '2 2' 
        hdr_out['WCSDIM'] = 2
        hdr_out['CTYPE1'] = 'LINEAR'
        hdr_out['CTYPE2'] = 'LINEAR'
        hdr_out['CRPIX1'] = 1. 
        hdr_out['CRVAL1'] = numpy.log10(newx[0])
        hdr_out['LTM1_1'] = 1. 
        hdr_out['LTM2_2'] = 1. 
        hdr_out['CDELT1'] = step
        #print newx[0], numpy.log10(newx[0]), step
        hdr_out['CD1_1'] = step
        hdr_out['CD2_2'] = 1.
        hdr_out['WAT0_001'] =  'system=equispec' 
        hdr_out['WAT1_001'] =  'wtype=linear label=Wavelength units=Angstroms'
        hdr_out['WAT2_001'] =  'wtype=linear'
        hdu_out.writeto(fname.replace(".fits", "1d.fits"), clobber=True)
      elif opts.fitstable:
        print newx,newy, newerr
        tbhdu = pyfits.BinTableHDU.from_columns([
          pyfits.Column(name='Arg', format='1D', disp='F8.3', array=newx),
          pyfits.Column(name='Fun', format='1D', disp='F8.3', array=newy),
          pyfits.Column(name='Var', format='1D', disp='F8.3', array=newerr)
          ])
        for i in range(1,1000):
          if hdr.__contains__("WAT2_%03d" % (i)):
            del hdr["WAT2_%03d" % (i)]
          else:
            break
        prihdu = pyfits.PrimaryHDU(header=hdr)
        tabhdr = tbhdu.header
        tabhdr['EXTNAME'] = 'DataVector'
        thdulist = pyfits.HDUList([prihdu, tbhdu])
        thdulist.writeto(fname.replace(".fits", "tab.fits"), clobber=True)
      '''
      # DEBUG: plot the spectrum
      import pylab
      pylab.subplot(2,1,1)
      for i in range(norder):
        pylab.plot(wd[i][t[i]],od_cont[i][t[i]])
      pylab.plot(newx,newy, lw=2, color="red")
      pylab.subplot(2,1,2)
      for i in range(norder):
        pylab.plot(wd[i][t2[i]],err_cont[i][t2[i]])
      pylab.plot(newx,newerr, lw=2, color="red")
      pylab.show()
      '''
    else:
        outfile = open(fname.replace(".fits", ".txt"), "w+")
        try:
          outfile.write("# OBJNAME = %s\n" % (hdr['OBJNAME']))
        except:
          outfile.write("# OBJNAME = %s\n" % (hdr['OBJECT']))
        try:
          outfile.write("# HJD = %s\n" % (hdr['HJD']))
        except:
          pass
        try:
          outfile.write("# VOBS = %s\n" % (hdr['VOBS']))
        except:
          pass
        try:
          outfile.write("# VORBIT = %s\n" % (hdr['VORBIT']))
        except:
          pass
        try:
          outfile.write("# VCORRECT = %s\n" % (hdr['VCORRECT']))
        except:
          outfile.write("# VCORRECT = %s\n" % (float(hdr['VHELIO'])-float(hdr['VOBS'])))
        for o in range(norder/nslices-1, -1, -1):
          for sl in range(nslices):
            ap = o*nslices+sl
            for n in range(len(wd[ap][t[ap]])):
              outfile.write("%f %f %f %d %d\n" % (wd[ap][t[ap]][n], od_cont[ap][t[ap]][n], sn[ap][t[ap]][n], o+1+offset, sl+1 ))
        outfile.close()

