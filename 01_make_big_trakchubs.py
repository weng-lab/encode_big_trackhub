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
from byAssayByBiosampleType import TrackhubDbByAssayByBiosampleType
from byAssayByFactor import TrackhubDbByAssayByFactor
from byCcREs import TrackhubDbByCcREs

sys.path.append(os.path.join(os.path.dirname(__file__), '../metadata/utils'))
from files_and_paths import Dirs, Urls, Datasets
from utils import Utils, eprint, AddPath, printt, printWroteNumLines
from metadataws import MetadataWS


class MegaTrackHub:
    def __init__(self, args, assembly, globalData):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData

        dataset = Datasets.byAssembly(assembly)
        self.mw = MetadataWS(dataset=dataset, host=Host)

    def run(self):
        self._makeHub()

        self.byBiosampleTypeOutput = ""
        if self.args.biosample:
            self.byBiosampleTypeOutput = TrackhubDbBiosampleType(self.args, self.assembly,
                                                                 self.globalData, self.mw).run()

        self.byAssayByBiosampleTypeOutput = ""
        if self.args.assay:
            self.byAssayByBiosampleTypeOutput = TrackhubDbByAssayByBiosampleType(self.args, self.assembly,
                                                                  self.globalData, self.mw).run()

        self.byAssayByFactorOutput = ""
        if self.args.factor:
            self.byAssayByFactorOutput = TrackhubDbByAssayByFactor(self.args, self.assembly,
                                                                   self.globalData, self.mw).run()

        self.byCcREsOutput = ""
        if self.args.ccREs:
            self.byCcREsOutput = TrackhubDbByCcREs(self.args, self.assembly,
                                                   self.globalData, self.mw).run()

        self.makeMainTrackDb()

    def makeMainTrackDb(self):
        fnp = os.path.join(BaseWwwDir, self.assembly, 'trackDb.txt')
        Utils.ensureDir(fnp)
        with open(fnp, 'w') as f:
            f.write(self.byAssayByBiosampleTypeOutput)
            f.write(self.byBiosampleTypeOutput)
            f.write(self.byCcREsOutput)
            f.write(self.byAssayByFactorOutput)
        printWroteNumLines(fnp)

    def _makeHub(self):
        fnp = os.path.join(BaseWwwDir, 'hub.txt')
        with open(fnp, 'w') as f:
            f.write("""
hub ENCODE
shortLabel ENCODE Trackhub Test5
longLabel ENCODE Trackhub Test5
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
    parser.add_argument("--assembly", type=str, default="")

    assay_parser = parser.add_mutually_exclusive_group(required=False)
    assay_parser.add_argument('--assay', dest='assay', action='store_true')
    assay_parser.add_argument('--no-assay', dest='assay', action='store_false')
    parser.set_defaults(assay=True)

    biosample_parser = parser.add_mutually_exclusive_group(required=False)
    biosample_parser.add_argument('--biosample', dest='biosample', action='store_true')
    biosample_parser.add_argument('--no-biosample', dest='biosample', action='store_false')
    parser.set_defaults(biosample=True)

    ccREs_parser = parser.add_mutually_exclusive_group(required=False)
    ccREs_parser.add_argument('--ccREs', dest='ccREs', action='store_true')
    ccREs_parser.add_argument('--no-ccREs', dest='ccREs', action='store_false')
    parser.set_defaults(ccREs=True)

    factor_parser = parser.add_mutually_exclusive_group(required=False)
    factor_parser.add_argument('--factor', dest='factor', action='store_true')
    factor_parser.add_argument('--no-factor', dest='factor', action='store_false')
    parser.set_defaults(factor=True)

    return parser.parse_args()


def main():
    args = parse_args()

    assemblies = ["hg19", "mm10"]
    for assembly in assemblies:
        printt("************************", assembly)

        if 0:
            printt("loading globalData from API...")
            globalDataUrl = "http://api.wenglab.org/screenv10_python/globalData/0/" + assembly
            ws = requests.get(globalDataUrl)
            globalData = ws.json()
        else:
            printt("loading globalData from disk...")
            fnp = os.path.join(os.path.dirname(__file__), "lists",
                               "globalData." + assembly + ".json")
            with open(fnp) as f:
                globalData = json.load(f)
        printt("done")

        tdb = MegaTrackHub(args, assembly, globalData)
        tdb.run()
    outputGenomes(assemblies)

if __name__ == '__main__':
    main()
