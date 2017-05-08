#!/usr/bin/env python
import os
import scipy as sp
import numpy as np
from scipy import signal
import multiprocessing
from multiprocessing import Pool
import shutil
import time
import tempfile
import commands
import ConfigParser
import commands
import logging
import sys
import optparse

from Lightcurve import AsciiLightcurve as lc
import util
from dataio import *

'''
Import Lightcurve, util and dataio from Chelsea's FFITools package.
'''

def cmd_parse():
    p = optparse.OptionParser()

    p.add_option('--cfgfile', '-c', default='example.cfg',
                 help='the filename of configuration file that '
                      'contains parameter settings(Optional)')
    p.add_option('--overwrite', '-r', default=False, action='store_true',
                 help='the flag to allow replacing existing files')
    p.add_option('--nprocs', '-n', default=3,
                 help='the number of processors to use when run the code with multiprocessing, default is 3, when set to 1, turn off multiprocessing')
    p.add_option('--logfile', default='example.log',
                 help='The destination of the log file (this is disabled when using multiple cores since multiprocessing does not take file handeler, output to stdout for now), - outputs the logs to stdout')
    p.add_option('--debug', action='store_true', help='output debug level message')

    options, arguments = p.parse_args()
    return options


def ls_parse(infile, name, default=None):
    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        var = p.get('LS', name)
    except IOError:
        raise IOError, 'Can not find the input configure file %s' % infile
    except ConfigParser.MissingSectionHeaderError:
        raise ConfigParser.MissingSectionHeadeError, 'Can not find Section LS in configure file %s\n' % infile
    except ConfigParser.NoOptionError:
        if default is None:
            raise
        else:
            logger.debug("Use default value for %s" % name)
            return default

    return var


class vt_ls(object):
    def __init__(self, cfgfile=''):
        logger.info("parsing cfg file")
        self.minp = float(ls_parse(cfgfile, 'minp', 0.5))
        self.subsample = float(ls_parse(cfgfile, 'subsample', 0.1))
        self.maxp = float(ls_parse(cfgfile, 'maxp', 30.))
        self.npeaks = int(ls_parse(cfgfile, 'npeaks', 3))
        self.nclip = int(ls_parse(cfgfile, 'nclip', 5))
        self.nclipiter = int(ls_parse(cfgfile, 'nclipiter', 1))

        self.length = 3
        return

    def __call__(self, lcfile, lsanal='', replace=False):
        lsanalfile = util.getext(lcfile.name, lsanal, '')
        lspath = os.path.dirname(lsanal)
        if lspath == "":
            lspath = "./"
        modelpath = lspath
        phasepath = lspath
        if os.path.exists(lsanalfile) and (not replace):
            logger.warning("File %s exists, do not allow to replace, set --rp to force overwrite", lsanal)
            return
        elif (not os.path.exists(lcfile.name)):
            logger.warning("Input file %s do not exists, do nothing", lcfile.name)
            return

        else:
            cmdline = "vartools -i %s -header -readformat 1 %d %d 3 -LS %f %f %f %d 1 %s whiten clip %d %d" % (
            lcfile.name, lcfile.cols['jd'], lcfile.cols['ltflc'], self.minp, self.maxp, self.subsample, self.npeaks, \
            lspath, self.nclip, self.nclipiter)
            logger.info("excuting cmdline: %s", cmdline)
            status, output = commands.getstatusoutput(cmdline)
            # print output
            header, result = output.split('\n')
            print result
            newheader = " ".join(header.split()[1:self.length + 1])
            newlines = ""
            for i in xrange(self.npeaks):
                newlines += "%d " % (i + 1)
                newlines += " ".join(result.split()[1 + i * self.length:(1 + self.length * (i + 1))]) + "\n"
            fout = open(lsanalfile, mode="w")
            fout.write("#")
            fout.write(newheader)
            fout.write("\n")
            fout.write(newlines.rstrip())
            fout.close()
        return


def lsmulti(args):
    lsengin = args[0]
    lcfile = args[1]
    lsanal = args[2]
    replace = args[3]
    lsengin(lcfile, dftanal=lsanal, replace=replace)
    return


def main(options, Multi=False):
    '''Documentation for runls module'''
    replace = options.overwrite
    nps = int(options.nprocs)
    lcfile = lc(cfgfile=options.cfgfile)

    infile = ls_parse(cfgfile, 'infile', '')
    if infile == '':
        logger.debug("Can not find input file from configure file, try to input from a input list from configure file")

    inlist = ls_parse(cfgfile, 'inlist', '')
    inpath = ls_parse(cfgfile, 'inpath', '')
    lsanal = ls_parse(cfgfile, 'analfile', '.lsanal')
    outpath = ls_parse(cfgfile, 'outpath', inpath)
    lsanal = outpath + lsanal
    colid = int(ls_parse(cfgfile, 'colid', 1))
    coljd = int(ls_parse(cfgfile, 'coljd', 1))
    colmag = int(ls_parse(cfgfile, 'colmag', 2))

    if (not infile == ''):
        lcfile.name = inpath + infile
        lcfile.cols['jd'] = coljd
        lcfile.cols['ltflc'] = colmag
        lsengin = vt_ls(cfgfile)
        lsengin(lcfile, lsanal=lsanal, replace=replace)
    else:
        names = []
        readcolumn(names, colid, inpath + inlist, datformat='str')
        lsengin = vt_ls(cfgfile)
        if (Multi):

            pool = Pool(processes=nps)
            lcfilelist = []
            for i in xrange(len(names)):
                lcfile = lc(cfgfile=options.cfgfile)
                lcfile.name = inpath + names[i]
                lcfile.cols['jd'] = coljd
                lcfile.cols['ltflc'] = colmag
                lcfilelist.append(lcfile)
            TASK = [([lsengin, lcfilelist[i], lsanal, replace]) for i in range(len(names))]
            pool.map(lsmulti, TASK)
        else:
            for x in names:
                # print x
                lcfile.cols['jd'] = coljd
                lcfile.cols['ltflc'] = colmag
                lcfile.name = inpath + x
                lsengin(lcfile, lsanal=lsanal, replace=replace)


if __name__ == '__main__':
    options = cmd_parse()
    if options.debug:
        loggerlevel = logging.DEBUG
    else:
        loggerlevel = logging.INFO

    cfgfile = options.cfgfile
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if int(options.nprocs) == 1:
        Multi = False
        logger = logging.getLogger(__name__)
        if options.logfile == '-':
            fh = logging.StreamHandler(stream=sys.stdout)
        else:
            fh = logging.FileHandler(options.logfile)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(logging.ERROR)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        logger.addHandler(fh)

    else:
        Multi = True
        logger = multiprocessing.get_logger()
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    logger.setLevel(loggerlevel)

    main(options, Multi=Multi)
