#!/usr/local/bin/python2.7
# encoding: utf-8
'''
digital_noise_image -- Script produces a digital noise image from DP an SP rounding differences

digital_noise_image is a main program using numarray and matplotlib

It defines classes_and_methods

@author:     andreas.wicenec@icrar.org

@copyright:  2015 ICRAR/UWA. All rights reserved.

@license:    Apache

@contact:    andreas.wicenec@icrar.org
@deffield    updated: Updated
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2015-09-15'
__updated__ = '2015-09-15'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def img_rounding_noise(psize=100):
    """
    Generate digital noise image from difference of SP and DP arrays.
    The SP array is converted to DP and then the difference is calculated.
    This effect is due to rounding differences between SP and DP.
    
    INPUTS:    
        psize:       int, size of image in pixels/side
    
    RETURNS:    
        the image as a numpy array (psize,psize).dtype=np.float64
    """
    na=np.array(np.random.sample(psize**2)*100,dtype=np.int8) # random integers [0..100[
    r32 = np.array(na, dtype=np.float32)/3 # produce SP array
    r64 = np.array(na, dtype=np.float64)/3 # produce DP array
    r32_64=np.array(r32,dtype=np.float64) # convert SP array to DP
    d64 = r64 - r32_64 # difference between the two (should be == 0!)
    d32 = np.array(d64,dtype=np.float32)
    img = d32.reshape(psize,psize) # make 2-D
    img += np.abs(img.min()) # shift to all positive
    return img


def main(argv=None): 
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by andreas.wicenec@icrar.org on %s.
  Copyright 2015 ICRAR/UWA. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-s", "--size", type=int, dest="psize", default=100, help="size of image pixels/side [default: %(default)s]")
        parser.add_argument("-v", "--verbose", dest="verbose", default=0, action="count", help="set verbosity level [default: %(default)s]")
        parser.add_argument("-m", "--method", dest="method", default='img_rounding_noise', help="set method to be used [default: %(default)s]")
        parser.add_argument("-l", "--list", dest="lm", action="count", default=False, help="list methods [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        # parser.add_argument(dest="paths", help="paths to folder(s) with source file(s) [default: %(default)s]", metavar="path", nargs='+')

        # Process arguments
        args = parser.parse_args()

        psize = args.psize
        verbose = args.verbose
        method = args.method
        
        if verbose > 0:
            print("Verbose mode on")
            if psize:
                print("Pixel size: %d" % psize)
            else:
                print("Default size 100")

        if args.lm:
            methods = filter(lambda x:x[0:4]=='img_',globals().keys())
            for m in methods:
                print m.__doc__
            return 0

        img_producer = globals()[method]
        img = img_producer(psize=psize)
        print "Maximum difference: %f" % np.max(img)
        plt.imshow(img)
        plt.show()
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-v")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'digital_noise_image_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())