#!/usr/bin/env python2

from __future__ import print_function

import sys
import json
import os
import re
import argparse
import requests
from collections import OrderedDict, defaultdict
from multiprocessing import Process, Value, Lock, Manager, Pool
from joblib import Parallel, delayed

from helpers.tracks import Tracks, Parent
import helpers.helpers as Helpers
from paths import Host, BaseWwwDir, BaseWwwTmpDir
from byBiosampleType import TrackhubDbBiosampleType
from byAssayByBiosampleType import TrackhubDbByAssayByBiosampleType
from byAssayByFactor import TrackhubDbByAssayByFactor
from byCcREs import TrackhubDbByCcREs
from byOrganSlim import TrackhubDbByOrganSlim

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from files_and_paths import Dirs, Urls, Datasets
from utils import Utils, eprint, AddPath, printt, printWroteNumLines
from metadataws import MetadataWS

class Counter(object):
    # https://github.com/davidheryanto/etc/blob/master/python-recipes/parallel-joblib-counter.py
    def __init__(self, manager, initval=0):
        self.val = manager.Value('i', initval)
        self.lock = manager.Lock()

    def increment(self, val = 1):
        return self.add_pre(val)

    # pre-fix add a given val
    def add_pre(self, val = 1):
        with self.lock:
            self.val.value += val
            return self.val.value

    # post-fix add a given val
    def add_post(self, val = 1):
        with self.lock:
            ret = self.val.value
            self.val.value += val
            return ret

    def value(self):
        with self.lock:
            return self.val.value


class MegaTrackHub:
    def __init__(self, args, assembly, globalData, priority):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData
        self.priority = priority

        dataset = Datasets.byAssembly(assembly)
        self.mw = MetadataWS(dataset=dataset, host=Host)

    def run(self):
        self._makeHub()

        args = {"args": self.args,
                "assembly": self.assembly,
                "globalData": self.globalData,
                "mw": self.mw,
                "priority": self.priority}

        def runner(arg, klass):
            if not arg:
                return ""
            return klass(**args).run()

        self.typs = [("ccREs", TrackhubDbByCcREs),
                     ("organSlim", TrackhubDbByOrganSlim),
                     ("factor", TrackhubDbByAssayByFactor),
                     ("assay", TrackhubDbByAssayByBiosampleType),
                     ("biosample", TrackhubDbBiosampleType)]

        self.out = {}
        for typ, f in self.typs:
            self.out[typ] = runner(getattr(self.args, typ), f)

        self.makeMainTrackDb()

    def makeMainTrackDb(self):
        fnp = os.path.join(BaseWwwDir, self.assembly, 'trackDb.txt')
        Utils.ensureDir(fnp)
        with open(fnp, 'w') as f:
            for typ, _ in self.typs:
                f.write(self.out[typ])
        printWroteNumLines(fnp)

    def _makeHub(self):
        fnp = os.path.join(BaseWwwDir, 'hub.txt')
        with open(fnp, 'w') as f:
            f.write("""
hub ENCODE
shortLabel ENCODE Trackhub Test6
longLabel ENCODE Trackhub Test6
genomesFile genomes.txt
email zhiping.weng@umassmed.edu
descriptionUrl http://encodeproject.org/
""")
        printWroteNumLines(fnp)

def outputGenomes(assemblies):
    fnp = os.path.join(BaseWwwDir, 'genomes.txt')

    defaultPoses = {"hg19": "chr12:121374959-121481905",
                    "mm10": "chr2:163423234-163655010"}

    with open(fnp, 'w') as f:
        for assembly in assemblies:
            defaultPos = defaultPoses.get(assembly, "")

            f.write("""
genome {assembly}
trackDb {assembly}/trackDb.txt
defaultPos {defaultPos}
""".format(assembly = assembly,
           defaultPos = defaultPos))
    printWroteNumLines(fnp)


def testHub():
    printt("checking hub...")
    cmds = ["/data/common/tools/ucsc.v350/hubCheck",
            "-noTracks",
            os.path.join(BaseWwwDir, 'hub.txt')]
    printt(Utils.runCmds(cmds))

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

    organSlim_parser = parser.add_mutually_exclusive_group(required=False)
    organSlim_parser.add_argument('--organSlim', dest='organSlim', action='store_true')
    organSlim_parser.add_argument('--no-organSlim', dest='organSlim', action='store_false')
    parser.set_defaults(organSlim=True)

    return parser.parse_args()


def main():
    args = parse_args()

    manager = Manager()
    result = manager.dict()
    priority = Counter(manager, 0)

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
            fnp = os.path.join(os.path.dirname(__file__), "../lists",
                               "globalData." + assembly + ".json")
            with open(fnp) as f:
                globalData = json.load(f)
        printt("done")

        tdb = MegaTrackHub(args, assembly, globalData, priority)
        tdb.run()
    outputGenomes(assemblies)
    testHub()

if __name__ == '__main__':
    main()
