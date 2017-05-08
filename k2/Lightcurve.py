#!/usr/bin/env p    @staticmethodython
import math
import numpy as np
from dataio import readcolumn
import os
import ConfigParser
import logging
class Lightcurve(object):

    def __init__(self, name=''):
        # Need to redesign this for different appertures.
        self.name = name
        self.data = {'jd': [], 'rlc': [], 'ltflc': [],
                     'x': [], 'y': [], 'bg': [], 'bgerr': [], 'cadence': []}
        self.labels = {'jd': 'BJD', 'rlc': 'Raw LC (Magnitude)',
                       'ltflc': 'DLC (Magnitude)', 'x': 'X (Pixel)',
                       'y': 'Y (Pixel)', 'bg': 'Background',
                       'bgerr': 'Background Err', 'cadence': 'Cadence'}
        self.logger = logging.getLogger(Lightcurve.__name__)
        self.logger.addHandler(logging.NullHandler())
        # self.ap = [True,True,True]
        return
    
    def __getstate__(self):
        d = self.__dict__.copy()
        if 'logger' in d.keys():
            d['logger'] = d['logger'].name
        return d

    def __setstate__(self, d):
        if 'logger' in d.keys():
            d['logger'] = logging.getLogger(d['logger'])
        self.__dict__.update(d)

    @staticmethod
    def create_lightcurve(infile,cfgfile='example.cfg', useridentifier=None):
        fin = open(infile,mode='r')
        identifier = fin.read(8)
        fin.close() 
        for c in [FITSLC, HDFLC, AsciiLightcurve]:
            if c.istype(identifier):
                return c(name=infile,cfgfile=cfgfile)
        if useridentifier:
            for c in [FITSLC, HDFLC, AsciiLightcurve]:
                if c.istype(useridentifier):
                    return c(name=infile,cfgfile=cfgfile)

        raise NotImplementedError, 'File format for input file %s is not recognized' % infile 

    def load_from_file(self, label='all'):
        # load data from a file, allow user to specify load all or specific columns.
        raise NotImplementedError
        return

    def write_to_file(self, outfile):
        # Output lightcurve to a file as designed.
        raise NotImplementedError
        return

    def gen_jd(self):
        self.check_data(label='cadence')
        self.data['jd'] = self.data['cadence'] / 48. + 2457827.0
        return

    def check_data(self, label='rlc'):
        if len(self.data[label]) == 0:
            try:
                self.load_from_file(label=label)
            except (IndexError, KeyError):
                raise IndexError, "Column %s is required, but cannot found in file %s" % (label, self.name)

        return 

    def detrend_lc(self, detrend_func):
        required_cols = detrend_func.required_cols
        for col in required_cols:
            self.check_data(col)
        self.data['ltflc'] = detrend_func.detrend(self.data)

        return

    def plot(self, show=True, label='rlc'):
        # only works for the sub
        import matplotlib
        from matplotlib import pyplot as plt
        self.check_data('jd')
        self.check_data(label)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if self.data['jd'][0] > 2454000:
            ax.plot(self.data['jd'] - 2454000, self.data[label], 'k.')
            ax.set_xlabel("BJD-2454000")
        else:
            ax.plot(self.data['jd'], self.data[label], 'k.')
            ax.set_xlabel("BJD")
        
        ax.set_ylabel(self.labels[label])
        yformatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(yformatter)
        ax.grid(True)
        if show:
            plt.show()
        else:
            pngfile = os.path.splitext(self.name)[0] + '.png'
            plt.savefig(pngfile)
        return


class AsciiLightcurve(Lightcurve):
    def __init__(self, cfgfile='example.cfg', name=''):
        super(AsciiLightcurve, self).__init__(name)
        self.logger = logging.getLogger(AsciiLightcurve.__name__)
        self.logger.addHandler(logging.NullHandler())
        p = ConfigParser.RawConfigParser()
        p.read(cfgfile)
    
        try: 
            self.cols = eval(p.get('LC', 'lcformat'))
            self.logger.debug("lcformat %s", str(self.cols))
        except ConfigParser.NoOptionError:
            self.cols = {'jd': 2, 'rlc': 8, 'ltflc': 12,
                        'x': 3, 'y': 4, "bg": 12, "bgerr": 12, 'cadence': 2}
            self.logger.warning("Use default LC format, type is %s", str(self.cols))  
        except ConfigParser.NoSectionError:
            raise ConfigParser.NoSectionError, "No Section LC in configure file %s" % (cfgfile)
        
        self.outext = ''
    
    @staticmethod 
    def istype(identifier):

        if all(ord(c) <128 for c in identifier): 
            return True
        return False

    def load_from_file(self, label='all'):
        # load data from a file, allow user to specify load all or specific columns.
        with open(self.name, mode='r') as fin:
            headerline = fin.readline()
        
        if headerline.startswith('#'):
            newdict = dict(zip(headerline.split()[1:],np.arange(len(headerline.split())-1)+1))
            self.cols.update(newdict)
            self.logger.debug("Find header line in the file, will try to use columns in the header instead, %s\n", str(self.cols))
        else: 
            self.logger.debug("File has no header line, use the default format or the format find in configure file, %s\n", str(self.cols))

        if label == 'all':
            readcolumn(self.data['jd'], int(self.cols['jd']), self.name)
            self.data['jd'] = np.array(self.data['jd'])
            readcolumn(self.data['rlc'], int(self.cols['rlc']), self.name)
            self.data['rlc'] = np.array(self.data['rlc'])
            try:
                readcolumn(self.data['ltflc'], int(self.cols['ltflc']), self.name)
                self.data['ltflc'] = np.array(self.data['ltflc'])
            except (IndexError, KeyError): 
                self.logger.warning("ltflc for %s not created yet, you might want to detrend first", self.name)
            try:
                readcolumn(self.data['x'], int(self.cols['x']), self.name)
                self.data['x'] = np.array(self.data['x'])
                readcolumn(self.data['y'], int(self.cols['y']), self.name)
                self.data['y'] = np.array(self.data['y'])
            except (IndexError, KeyError): 
                self.logger.warning("coordinate columns x, " \
                      "y for %s not created yet, you might want " \
                      "to detrend first", self.name)
            try:
                readcolumn(self.data['bg'], int(self.cols['bg']), self.name)
                self.data['bg'] = np.array(self.data['bg'])
                readcolumn(self.data['bgerr'], int(self.cols['bgerr']), self.name)
                self.data['bgerr'] = np.array(self.data['bgerr'])
            except (IndexError, KeyError): 
                self.logger.warning("Warning: background columns bg, bgerr for %s not " \
                      "created yet, you might want to detrend first", self.name)
        else:
        
            readcolumn(self.data[label], int(self.cols[label]), self.name)
            self.data[label] = np.array(self.data[label])


    def write_to_file(self, outfile, replace=True, outheader="# jd rlc"):
        # Output lightcurve to a file as designed.
        if os.path.exists(outfile) and (not replace):
            self.logger.warning("outfile %s is already there, not allowed to replace using current cofiguration", outfile)
            return
        self.logger.info("outfile header is %s", outheader)
        with open(outfile, mode='w') as fout:
            fout.write(outheader)
            fout.write('\n')	    
            for i in xrange(len(self.data['jd'])):
                if np.isnan(self.data['rlc'][i]):
                    continue
                for label in outheader.split()[1:]:
                    fout.write("%15.7f " % (self.data[label][i]))
                fout.write("\n") 

class HDFLC(Lightcurve):
    # TBD. need to merge from old branch
    IDENTIFIER = "\211HDF\r\n\032\n"
    def __init__(self,name=''):
        super(HDFLC, self).__init__(name)
        self.logger = logging.getLogger(HDFLC.__name__)
        self.logger.addHandler(logging.NullHandler())
        return
    
    @staticmethod
    def istype(identifier):
        if identifier == HDFLC.IDENTIFIER:
            return True
        return  False
    
    def load_from_file(self, label='all'):
        raise NotImplementedError

    def write_to_file(self, outfile, replace=True, outheader="# jd rlc"):
        raise NotImplementedError

class FITSLC(Lightcurve):
    IDENTIFIER = "SIMPLE  "
    BJD = 2454000
    #Z0 = 24.2415 #Kepler
    Z0 = 20.17 #tess
    def __init__(self, name='', cfgfile='example.cfg'):
        super(FITSLC, self).__init__(name)
        self.logger = logging.getLogger(FITSLC.__name__)
        self.logger.addHandler(logging.NullHandler())

    @staticmethod
    def istype(identifier):
        if identifier == FITSLC.IDENTIFIER:
            return True
        return False
    def load_from_file(self, label='all'):
        # load data from a file, allow user to specify load all or specific columns.
        import astropy
        import astropy.io.fits as fits 
        with fits.open(self.name) as hdulist:
            tbdata = hdulist[1].data
            
            BJDREFI = hdulist[1].header['BJDREFI']
            BJDREFF = hdulist[1].header['BJDREFF']
     
            self.data['jd'] = np.array(tbdata.field('Time')) + BJDREFI + BJDREFF - FITSLC.BJD
            self.logger.debug("lenth of data is %d", len(self.data['jd']))
            self.data['rlc'] = FITSLC.Z0 - 2.5*np.log10(np.array(tbdata.field('PDCSAP_FLUX')))
            try:
                TESSmag=hdulist[0].header['TESSMAG']
                self.data['rlc']-=(np.nanmedian(self.data['rlc'])-TESSmag)
            except KeyError:
                pass
            #self.data['rlc'] = FITSLC.Z0 - 2.5*np.log10(np.array(tbdata.field('SAP_FLUX')))
            self.data['x'] = np.array(tbdata.field('MOM_CENTR1'))
            self.data['y'] = np.array(tbdata.field('MOM_CENTR2')) 
            self.data['cadence'] = np.array(tbdata.field('CADENCENO')) 
            # ?self.data['bg'] = tbdata.field['MOM_CENTR2'] 
            # ?self.data['bgerr'] = tbdata.field['MOM_CENTR2'] 

    def write_to_file(self, outfile, replace=True, outheader="# jd rlc"):
        # Output lightcurve to a file as designed.
        if os.path.exists(outfile) and (not replace):
            print "outfile %s is already there, not allowed to replace using current cofiguration" % outfile
            return
        self.logger.info("outfile header is %s", outheader)
        with open(outfile, mode='w') as fout:
            fout.write(outheader)
            fout.write('\n')	    
            for i in xrange(len(self.data['jd'])):
                if np.isnan(self.data['rlc'][i]):
                    continue
                for label in outheader.split()[1:]:
                    fout.write("%15.7f " % (self.data[label][i]))
                fout.write("\n") 

