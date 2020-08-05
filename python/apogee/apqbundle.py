#!/usr/bin/env python

"""APQBUNDLE.PY - This bundles the raw APOGEE frames and optionally compresses them.

The data will be output in this format:
 HDU0: header but NO data
 HDU1: header, read1 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
 HDU2: header, read2 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
 and so on for all the reads.

If /compress is set then APZIP will be used to compress the raw
APOGEE data.  Use APUNZIP to uncompress the data.
"""

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@noao.edu>'
__version__ = '20180922'  # yyyymmdd                                                                                                                           

import os
import numpy as np
import warnings
from astropy.io import fits
from astropy.table import Table, Column
from glob import glob
from apogee.utils import yanny, apload
from sdss_access.path import path

# Ignore these warnings, it's a bug
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def apbundle(framenum,indir=None,outdir=None,compress=False,clobber=True,outfiles=None):
    """This bundles the raw APOGEE frames and optionally compresses them.

    The data will be output in this format:
     HDU0: header but NO data
     HDU1: header, read1 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
     HDU2: header, read2 image  as UNSIGNED INTEGERS (BITPIX=16 or UINT)
     and so on for all the reads.

    If /compress is set then APZIP will be used to compress the raw
    APOGEE data.  Use APUNZIP to uncompress the data.

    Parameters
    ----------
    frameid : int or string
            The 8 character frame ID.
    indir : str, optional
          The input directory with the apRaw files.  If not given then
            the data model will be used to figure out the directory.
    outdir : str, optional
            The output directory for the apR-[abc]-frameid.fits files.
            If not given then the data model will be used.s
    compress : bool, optional
             Compress the output files.  False by default.
    clobber : bool, optional
            Overwrite files if they already exist. True by defalt.


    Returns
    -------
    outfiles : list
            The output files names.  Normally this is the bundles
            files, but if /compress is set then this is the
            compressed filenames.
            
    The raw bundles APOGEE files are names outdir/apR-[abc]-frameid.fits
    If /compress is set then compressed files will also be created and
    have names of outdir/apR-[abc]-frameid.apz.

    Example
    -------

    .. code-block:: python

        apbundle('00000117')

    """

    npix = 2048

    # ID8
    frameid = ("%08d" % np.int(framenum))
    mjd = apload.cmjd(np.int(framenum))

    # Get apRaw directory
    sdss_path = path.Path()
    rawdir = sdss_path.dir('apRaw',num=framenum,read=1,mjd=mjd)
    if observatory=='apo':
        prefix = 'ap'
    else:
        prefix = 'as'
    archive_dir = sdss_path.dir(prefix+'R',mjd=mjd,chip='a',num=1)
    if outdir is None: outdir = archive_dir
    
    # Find the files
    files = np.sort(glob(rawdir+'/'+prefix+'Raw-'+frameid+'-???.fits'))
    nfiles = len(files)
    if nfiles==0:
        raise Exception("No files for "+frameid)
    num = np.zeros(nfiles,dtype=int)
    for i in range(nfiles):
        base = os.path.basename(files[i])
        # apRaw-28190009-059.fits
        num[i] = np.int(base[15:18])
    maxind = np.argmax(num)
    print("Nfiles = "+str(nfiles)+" for "+frameid)

    # Read in the FIRST frame and initialize the output files
    #  this way the timestamp will be from the BEGINNING of
    #  the exposure
    minind = np.argmin(num)
    firstfile = files[minind]
    #FITS_READ,firstfile,im,head,exten=0,/no_abort,message=message
    im,head = fits.getdata(firstfile,0,header=True)

    chiptag = ['a','b','c']
    outfiles = []
    for ch in chiptag:
        o = os.path.join(outdir,os.path.basename(sdss_path.full(prefix+'R',mjd=mjd,chip=ch,num=frameid).replace('apz','fits')))
        outfiles.append(o)

    # Final Xsize
    shape = im.shape
    xsizemap = {6144: 2048, 8192: 2048+512}
    if shape[1] not in xsizemap:
        raise Exception('Non-standard Xsize = '+str(sz[1])+' not supported')
    xsize = xsizemap[shape[1]]

    # Check if the output files exist
    test = [os.path.exists(o) for o in outfiles]
    if (np.sum(test)>0) & (clobber is False):
      raise Exception('Output files exist already.  Set /clobber to overwrite')

    print('Writing bundled file to '+" ".join(outfiles))

    # Initialize the three bundled files
    ahead = head.copy()
    ahead["SIMPLE"] = "T"
    ahead["BITPIX"] = 16
    ahead["NAXIS"] = 0
    ahead["CHIP"] = "a"
    ahead["NREAD"] = nfiles
    remkeys = ["NAXIS1","NAXIS2","PCOUNT","GCOUNT","CHECKSUM","DATASUM","BZERO","BSCALE"]
    for k in remkeys:
        if ahead.get(k) is not None: ahead.remove(k)
    fits.writeto(outfiles[0],None,ahead,checksum=True,overwrite=True)

    bhead = ahead.copy()
    bhead["CHIP"] = "b"
    remkeys = ["CHEKSUM","DATASUM"]
    for k in remkeys:
        if bhead.get(k) is not None: bhead.remove(k)
    fits.writeto(outfiles[1],None,bhead,checksum=True,overwrite=True)

    chead = ahead
    chead["CHIP"] = "c"
    remkeys = ["CHEKSUM","DATASUM"]
    for k in remkeys:
        if chead.get(k) is not None: chead.remove(k)    
    fits.writeto(outfiles[2],None,chead,checksum=True,overwrite=True)

    # Read in the frames
    for fpos in range(nfiles):
        # to use the fits_test_checksum, the data has to be read as a signed int first
        # this is because the pyfits checksum is done on signed values
        im,head = fits.getdata(files[fpos],0,header=True)
        #FITS_READ,files[fpos],im,head,exten=0,/no_abort,/no_unsigned,/noscale,/no_pdu,message=message
        #res = FITS_TEST_CHECKSUM(head, im, errmsg=errmsg)
        #if res eq -1 then begin
        #  # the checksum doesn't match -> send a warning and keep going
        #  message,/con,'BAD checksum for file '+files[fpos]
        #  return

        # Convert back to uint
        bzero = head.get("ZERO")
        if bzero is None: bzero=0
        if bzero!=0:
            im = np.uint(im)-np.uint(bzero[0])
        bscale = head.get("BSCALE")
        if bscale is None: bscale=1
        if bscale!=1:
            im /= bscale[0]

        # The image is: red, green, blue
        # if the 5th reference output is also there then
        # this is: green-window (not used), blue-ref, green-ref, red-ref

        # The RAW images:
        # Basic orientation is long wavelengths on the left, to shorter to the
        # right (higher column number).  readout direction is:
        #         Red                       Green                 Blue           WG RefB RefG RefR
        #|---->|<----|--->|<----||---->|<----|--->|<----||---->|<----|--->|<----||---->|<----|--->|<----|


        # Chip a, RED
        #--------------
        im1 = im[:,0:npix]
        # Add the 5th output, 4th column
        #  need to flip in X-direction so all three chips
        #  have the reference reading out left-to-right (--->)
        if xsize==2560:
            refim1 = im[:,6144+3*512:6144+4*512]
            refim1 = np.flip(refim1,axis=1)  # flip X-direction
            im1 = np.hstack((im1,refim1))
        ahead = head.copy()
        ahead.set("XTENSION","IMAGE",before="SIMPLE")
        bitpix = ahead.get("BITPIX")
        if bitpix is not None: ahead["BITPIX"] = bitpix
        ahead["NAXIS"] = 2
        ahead.set("NAXIS1",xsize,after="NAXIS")
        ahead.set("NAXIS2",len(im1[0,:]),after="NAXIS1")
        ahead.set("PCOUNT",0,after="NAXIS2")
        ahead.set("GCOUNT",1,after="PCOUNT")
        # delete SIMPLE if present, only allowed in PDU
        remkeys = ["SIMPLE","CHECKSUM","DATASUM"]
        for k in remkeys:
            if ahead.get(k) is not None: ahead.remove(k)
        ahead["CHIP"] = "a"
        ahead["NREAD"] = nfiles
        fits.append(outfiles[0],im1,ahead,checksum=True)

        # re-write the bzero and bscale to add the '/' that mwrfits automatically adds
        # so as to not mess up the checksum
        #if bzcount>0: sxaddpar,ahead,'BZERO',bzero
        #if bscount>0: sxaddpar,ahead,'BSCALE', bscale

        # Chip b, GREEN
        #---------------
        im2 = im[:,npix:2*npix-1]
        # Add the 5th output, 3rd column
        #  the reference readout direction is already left-to-right
        if xsize==2560:
            refim2 = im[:,6144+2*512:6144+3*512]
            im2 = np.hstack((im2,refim2))
        bhead = ahead.copy()
        bhead["NAXIS"] = 2
        bhead.set("NAXIS1",xsize,after="NAXIS")
        bhead.set("NAXIS2",len(im2[0,:]),after="NAXIS1")
        bhead.set("PCOUNT",0,after="NAXIS2")
        bhead.set("GCOUNT",1,after="PCOUNT")
        fits.append(outfiles[1],im2,bhead,checksum=True)

        # Chip c, BLUE
        #--------------
        im3 = im[:,2*npix:3*npix-1]
        # Add the 5th output, 2nd column
        #  need to flip in X-direction so all three chips
        #  have the reference reading out left-to-right (--->)
        if xsize==2560:
            refim3 = im[:,6144+512:6144+2*512]
            refim3 = np.flip(refim3,axis=1)  # flip X-direction
            im3 = np.hstack((im3,refim3))
        chead = ahead.copy()
        chead["NAXIS"] = 2
        chead.set("NAXIS1",xsize,after="NAXIS")
        chead.set("NAXIS2",len(im3[0,:]),after="NAXIS1")
        chead.set("PCOUNT",0,after="NAXIS2")
        chead.set("GCOUNT",1,after="PCOUNT")
        fits.append(outfiles[2],im3,chead,checksum=True)


        # Compress
        #----------
        if compress is True:
            print,'Compressing ',addslash(outdir)+file_basename(apogee_filename('R',num=frameid,chip=chiptag),'.apz')+'.fits'
            apzip(outfiles)

            # Return compressed filenames
            bundle_files = outfiles
            outfiles = file_dirname(bundle_files)+'/'+file_basename(bundle_files,'.fits')+'.apz'

    dt = systime(1)-t0
    print('dt = '+str(dt)+' sec')

    # 17 sec for a 60 read file.
    # 204 sec = 3.4 min for 60 read file WITH compression
