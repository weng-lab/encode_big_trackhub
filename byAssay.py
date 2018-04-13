from __future__ import print_function

import sys
import json
import os
import re
import argparse
import requests
from itertools import groupby
from collections import OrderedDict, defaultdict
from joblib import Parallel, delayed
import StringIO

from helpers.tracks import Tracks, Parent, LookupActive
import helpers.helpers as Helpers
from paths import Host, BaseWwwDir, BaseWwwTmpDir

sys.path.append(os.path.join(os.path.dirname(__file__), '../metadata/utils'))
from files_and_paths import Dirs
from utils import Utils, eprint, AddPath, printt, printWroteNumLines
from metadataws import MetadataWS
from files_and_paths import Urls

# from http://stackoverflow.com/a/19861595
import copy_reg
import types

def _reduce_method(meth):
    return (getattr, (meth.__self__, meth.__func__.__name__))
copy_reg.pickle(types.MethodType, _reduce_method)

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

class TrackhubDbByAssay:
    def __init__(self, args, assembly, globalData, mw):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData
        self.mw = mw

        self.expsByAssay= [("Chromatin Accessibility", "chromatin_accessibility",
                            self.mw.dnases_useful),
                           ("Histone Modifications", "histone_modifications",
                            self.mw.chipseq_histones_useful),
                           ("Transcription", "transcription",
                            self.mw.transcription_useful),
                           ("Transcription Factors", "transcription_factors",
                            self.mw.chipseq_tfs_useful)
        ]


        # assay x biosamepleType x biosamplesView

        self.byAssayBiosampleType = defaultdict(lambda: defaultdict(dict))
        self.subGroups = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.btToNormal = {}
        self.lookupByExp = {}

    def run(self):
        for title, assayAbbr, expsF in self.expsByAssay:
            exps = expsF()
            self._build(title, assayAbbr, exps)
        return self._makeMainTrackDb()

    def _build(self, assay_term_name, atn, exps):
        printt("building", assay_term_name, "...")
        def sorter(exp):
            return (exp.biosample_type)
        exps.sort(key = sorter)

        self.btToNormal[atn] = assay_term_name

        for biosample_type, exps in groupby(exps, sorter):
            exps = list(exps)
            bt = Helpers.sanitize(biosample_type)
            self.btToNormal[bt] = biosample_type

            fnp = os.path.join(BaseWwwTmpDir, self.assembly, "subtracks", atn, bt +'.txt')
            self.byAssayBiosampleType[atn][bt] = {
                "assay_term_name": assay_term_name,
                "atn": atn,
                "biosample_type": biosample_type,
                "bt": bt,
                "fnp": fnp,
                "exps": exps,
                "assembly": self.assembly
            }

        printt("making tracks and subtracks...")
        self._makeSubTracks()

    def _makeSubTracks(self):
        jobs = []
        for atn, btAndInfo in self.byAssayBiosampleType.iteritems():
            for bt, info in btAndInfo.iteritems():
                jobs.append(merge_two_dicts(info,
                                            {"lookupByExp": self.lookupByExp,
                                             "idx": len(jobs) + 1,
                                             "total": len(self.byAssayBiosampleType)}))

        Parallel(n_jobs=self.args.j)(delayed(outputAllTracksByBiosampleType)(job)
                                     for job in jobs)

    def _makeMainTrackDb(self):
        mainTrackDb = []

        priority = 0
        for atn, btAndInfo in self.byAssayBiosampleType.iteritems():
            priority += 1
            totalExperiments = sum([len(info["exps"]) for info in btAndInfo.values()])
            shortLabel = self.btToNormal[atn]
            longLabel = self.btToNormal[atn] + " (%s experiments)" % totalExperiments
            mainTrackDb.append("""
track super_{atn}
superTrack on show
priority {priority}
shortLabel {shortL}
longLabel {longL}
""".format(atn = atn,
           priority = priority,
           shortL=Helpers.makeShortLabel(shortLabel),
           longL=Helpers.makeLongLabel(longLabel)))

        outF = StringIO.StringIO()
        outF.write('\n'.join(mainTrackDb))
        for atn, btAndInfo in self.byAssayBiosampleType.iteritems():
            for bt, info in btAndInfo.iteritems():
                with open(info["fnp"]) as f:
                    outF.write(f.read())
                    outF.write('\n')
        return outF.getvalue()

def outputAllTracksByBiosampleType(info):
    subGroups = outputSubTrack(**info)
    info["subGroups"] = subGroups
    outputCompositeTrackByBiosampleType(**info)

def outputCompositeTrackByBiosampleType(assembly, assay_term_name, atn, biosample_type, bt,
                                        exps, fnp, idx, total, subGroups, lookupByExp):
    if not os.path.exists(fnp):
        raise Exception("missing " + fnp)

    subGroupsDict = {}
    for k in Helpers.SubGroupKeys:
        subGroupsDict[k] = {a[0]:a[1] for a in subGroups[k]}
    longLabel = assay_term_name + '_' + biosample_type + " (%s experiments)" % len(exps)

    print(subGroupsDict["assay"])

    if "chromatin_accessibility" == atn:
        subGroup1key = "biosample"
        subGroup2key = "age"
    elif "histone_modifications" == atn:
        subGroup1key = "biosample"
        subGroup2key = "label"
    elif "transcription_factors" == atn:
        subGroup1key = "biosample"
        subGroup2key = "label"
    elif "transcription" == atn:
        subGroup1key = "biosample"
        subGroup2key = "assay"
    else:
        subGroup1key = "donor"
        subGroup2key = "age"
    subGroup3key = "view"
    subGroup1 = Helpers.unrollEquals(subGroupsDict[subGroup1key])
    subGroup2 = Helpers.unrollEquals(subGroupsDict[subGroup2key])
    subGroup3 = Helpers.unrollEquals(subGroupsDict[subGroup3key])

    with open(fnp) as f:
        subtracks = f.read()

    actives = []
    # for expID in expIDs:
    #     if expID in lookupByExp:
    #         actives.append(lookupByExp[expID].isActive())
    isActive = any(t for t in actives)
    # if isActive:
    #     print("active biosample (composite):", btn)

    with open(fnp, 'w') as f:
        f.write("""
track {atn}_{bt}
parent super_{atn}
compositeTrack on
""".format(atn=atn,
           bt=bt))
        if isActive:
            f.write("visibility full\n")
        f.write("""shortLabel {shortL}
longLabel {longL}
type bigWig 9 +
maxHeightPixels 64:12:8
autoScale on
subGroup1 {subGroup1key} {subGroup1key} {subGroup1}
subGroup2 {subGroup2key} {subGroup2key} {subGroup2}
subGroup3 {subGroup3key} {subGroup3key} {subGroup3}
sortOrder {subGroup1key}=+ {subGroup2key}=+ {subGroup3key}=+
dimensions dimX={subGroup2key} dimY={subGroup1key}
dragAndDrop subTracks
hoverMetadata on
darkerLabels on
""".format(shortL=Helpers.makeShortLabel(biosample_type),
           longL=Helpers.makeLongLabel(longLabel),
           subGroup1key=subGroup1key,
           subGroup1=subGroup1,
           subGroup2key=subGroup2key,
           subGroup2=subGroup2,
           subGroup3key=subGroup3key,
           subGroup3=subGroup3
))
        f.write('\n' + subtracks)

    printWroteNumLines(fnp, idx, 'of', total)

def outputSubTrack(assembly, assay_term_name, atn, biosample_type, bt,
                   exps, fnp, idx, total, lookupByExp):
    actives = []
    # for expID in expIDs:
    #     if expID in lookupByExp:
    #         actives.append(lookupByExp[expID].isActive())
    isActive = any(t for t in actives)
    if isActive:
        print("active biosample:", btn)

    parent = Parent(atn + '_' + bt, isActive)

    tracks = Tracks(assembly, parent, (1 + idx) * 1000)
    for exp in exps:
        active = False
        expID = exp.encodeID
        cREs = {}
        if expID in lookupByExp:
            active = lookupByExp[expID].isActive()
            cREs = lookupByExp[expID].cREs
        tracks.addExp(exp, active, cREs)

    Utils.ensureDir(fnp)
    with open(fnp, 'w') as f:
        for line in tracks.lines():
            f.write(line)
    printWroteNumLines(fnp, idx, 'of', total)
    return tracks.subgroups()
