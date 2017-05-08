import ConfigParser
import optparse
import logging
import sys
import os
import multiprocessing
from multiprocessing import Pool

from Lightcurve import Lightcurve, AsciiLightcurve, FITSLC, HDFLC 
from Detrend import DetrendFunc, COSDetrend, EPDDetrend, FocusDetrend, PolyDetrend

def lc_parse(infile, name, default=None, **kwargs):
    if 'logger' in kwargs:
        logger = kwargs['logger']
    else: 
        logger = logging.getLogger(lc_parse.__name__)
    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        var = p.get('LC', name)
    except IOError:
        raise IOError, 'Can not find the input configure file %s' % infile
    except ConfigParser.MissingSectionHeaderError:
        raise ConfigParser.MissingSectionHeadeError, 'Can not find Section LC in configure file %s\n' % infile 
    except ConfigParser.NoOptionError:
        if default is None:
            raise
        else:
            logger.debug("Use default value for %s" % name)
            return default
    return var



def cmd_parse():
    p = optparse.OptionParser()
    p.add_option('--logfile', default='example.log', 
                 help='The destination of the log file, - outputs the logs to stdout')
    p.add_option('--overwrite', '-r', default=False, action='store_true',
                 help='the flag to allow replacing existing files')
    p.add_option('--nprocs', '-n', default=3, help='the number of processors to use when run the code with multiprocessing, default is 3, when set to 1, turn off multiprocessing')

    p.add_option('--cfgfile', '-c', default='example.cfg',
                 help='the configure file that the pipeline '
                      'will read the settings for fitsh from, '
                      'default is example.cfg')
    p.add_option('--debug', action='store_true', help='output debug level message')
    options, arguments = p.parse_args()
    return options


def detrendmulti(args):
    print args
    infile = args[0]
    outfile = args[1]
    cfgfile = args[2]
    replace = args[3]
    detrend_func = args[4]
    lc = Lightcurve().create_lightcurve(infile,cfgfile=cfgfile)
    logger.info("Detrending %s...", lc.name)
    jdflag = True
    if not jdflag:
        lc.gen_jd()
    lc.detrend_lc(detrend_func)
    logger.info("output to %s...", outfile)
    lc.write_to_file(outfile, replace=replace, outheader="# jd rlc ltflc") 
    return

def main(options, cfgfile="example.cfg", Multi=False):
    replace = options.overwrite
    nps = int(options.nprocs)
    method = lc_parse(cfgfile, 'method', 'cos')
    detrend_func = DetrendFunc.create(method, cfgfile=cfgfile) 
    jdflag = lc_parse(cfgfile, 'jdflag', True)
    infile = lc_parse(cfgfile, 'infile', '')
    if infile=='':
        logger.warning('can not find a input file in cfg file %s Section LC, will try to read inlist as the input', cfgfile) 
        inlist = lc_parse(cfgfile, 'inlist', None)
        
    inpath = lc_parse(cfgfile, 'inpath', '')
    outpath = lc_parse(cfgfile, 'outpath', inpath)
    outext = lc_parse(cfgfile, 'outext', '_ltf.lc')
    #logger.info('Inpath set to be %s, outpath set to be %s, infile set to be %s, inlist set to be %s', inpath, outpath, infile, inlist)    
    if not infile == '':
        lc = Lightcurve().create_lightcurve(inpath+infile,cfgfile=cfgfile) 
        if not jdflag:
            lc.gen_jd()
   
        lc.detrend_lc(detrend_func)
        outfile = outpath+os.path.basename(os.path.splitext(infile)[0])+outext
        lc.write_to_file(outfile, replace=replace, outheader="# jd rlc ltflc")
    elif not inlist == '':
        if Multi:
            pool = Pool(processes=nps)
            infilelist = []
            outfilelist = []
            for line in open(inpath+inlist).readlines():
                name = line.rstrip()
                infile = inpath+name 
                
                outfile = outpath+os.path.basename(os.path.splitext(name)[0]) + outext
                #print infile, outfile
                infilelist.append(infile)
                outfilelist.append(outfile)
            TASK=[([infilelist[i], outfilelist[i], cfgfile, replace, detrend_func]) for i in range(len(infilelist))]
            pool.map(detrendmulti, TASK)
        else:
            for line in open(inpath+inlist).readlines():
                name = line.rstrip()
                lc = Lightcurve().create_lightcurve(inpath+name,cfgfile=cfgfile) 
                logger.info("Detrending %s...", lc.name)
                if not jdflag:
                    lc.gen_jd()
   
                lc.detrend_lc(detrend_func)
                outfile = outpath+os.path.basename(os.path.splitext(name)[0]) + outext
                logger.info("output to %s...", outfile)
                lc.write_to_file(outfile, replace=replace, outheader="# jd rlc ltflc")
    return


if __name__ == '__main__':
    options = cmd_parse()
    if options.debug:
        loggerlevel = logging.DEBUG
    else:
        loggerlevel = logging.INFO
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if int(options.nprocs) == 1:
        Multi = False
        logger = logging.getLogger(__name__)
        if options.logfile == '-':
            fh = logging.StreamHandler(stream=sys.stdout)
        else:
            fh = logging.FileHandler(options.logfile)
        fh.setLevel(logging.DEBUG)
        # create console handler with a higher log level
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(logging.ERROR)
        # create formatter and add it to the handlers
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)
        logger.addHandler(ch)
    else:
        Multi = True
        logger = multiprocessing.get_logger()
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    logger.setLevel(loggerlevel)
    main(options, cfgfile=options.cfgfile, Multi=Multi)
