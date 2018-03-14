#!/usr/bin/env python2

from __future__ import print_function

import sys
import json
import os
import re
import argparse
import requests
from collections import OrderedDict, defaultdict
from joblib import Parallel, delayed
import StringIO

from helpers.tracks import Tracks, Parent
import helpers.helpers as Helpers
from paths import Host, BaseWwwDir, BaseWwwTmpDir
from byBiosampleType import TrackhubDbBiosampleType

sys.path.append(os.path.join(os.path.dirname(__file__), '../metadata/utils'))
from files_and_paths import Dirs
from utils import Utils, eprint, AddPath, printt, printWroteNumLines
from metadataws import MetadataWS
from files_and_paths import Urls


class MegaTrackHub:
    def __init__(self, args, assembly, globalData):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData

        self.mw = MetadataWS(host=Host)

    def run(self):
        self._makeHub()

        self.byBiosampleType = TrackhubDbBiosampleType(self.args, self.assembly, self.globalData)
        self.byBiosampleTypeOutput = self.byBiosampleType.run()
        
        self.makeMainTrackDb()
        
    def makeMainTrackDb(self):
        fnp = os.path.join(BaseWwwDir, self.assembly, 'trackDb.txt')
        Utils.ensureDir(fnp)
        with open(fnp, 'w') as f:
            f.write(self.byBiosampleTypeOutput)
        printWroteNumLines(fnp)
        
    def _makeHub(self):
        fnp = os.path.join(BaseWwwDir, 'hub.txt')
        with open(fnp, 'w') as f:
            f.write("""
hub ENCODE
shortLabel ENCODE Trackhub Test3
longLabel ENCODE Trackhub Test3
genomesFile genomes.txt
email zhiping.weng@umassmed.edu
descriptionUrl http://encodeproject.org/
""")
        printWroteNumLines(fnp)

def outputGenomes(assemblies):
    fnp = os.path.join(BaseWwwDir, 'genomes.txt')
    with open(fnp, 'w') as f:
        for assembly in assemblies:
            f.write("""
genome {assembly}
trackDb {assembly}/trackDb.txt
""".format(assembly = assembly))
    printWroteNumLines(fnp)
  


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-j', type=int, default=4)
    parser.add_argument("--assembly", type=str, default="hg19")
    return parser.parse_args()


def main():
    args = parse_args()

    assemblies = ["hg19", "mm10"]
    for assembly in assemblies:
        printt("************************", assembly)

        printt("loading globalData from API...")
        globalDataUrl = "http://api.wenglab.org/screenv10_python/globalData/0/" + assembly
        ws = requests.get(globalDataUrl)
        globalData = ws.json()
        printt("done")

        tdb = MegaTrackHub(args, assembly, globalData)
        tdb.run()
    outputGenomes(assemblies)

if __name__ == '__main__':
    main()
