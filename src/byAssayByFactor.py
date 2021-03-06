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

class TrackhubDbByAssayByFactor:
    def __init__(self, args, assembly, globalData, mw, priority):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData
        self.mw = mw
        self.priority = priority

        self.expsByAssay= [("TFs by Factor",
                            "tf_factors",
                            self.mw.chipseq_tfs_useful),
                           ("Histone by Mark",
                            "hm_by_marks",
                            self.mw.chipseq_histones_useful)
        ]

        # assay x factor x biosamplesView

        self.byAssayBiosampleType = defaultdict(lambda: defaultdict(dict))
        self.subGroups = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.labelNToNormal = {}
        self.lookupByExp = {}

    def run(self):
        for title, assayAbbr, expsF in self.expsByAssay:
            exps = expsF()
            self._build(title, assayAbbr, exps)
        return self._makeMainTrackDb()

    def _build(self, assay_term_name, atn, exps):
        printt("building", assay_term_name, "...")
        def sorter(exp):
            return (exp.label)
        exps.sort(key = sorter)

        self.labelNToNormal[atn] = assay_term_name

        for label, exps in groupby(exps, sorter):
            exps = list(exps)
            labelN = Helpers.sanitize(label)
            self.labelNToNormal[labelN] = label

            fnp = os.path.join(BaseWwwTmpDir, self.assembly, "subtracks", atn, labelN +'.txt')
            self.byAssayBiosampleType[atn][labelN] = {
                "assay_term_name": assay_term_name,
                "atn": atn,
                "label": label,
                "labelN": labelN,
                "fnp": fnp,
                "exps": exps,
                "assembly": self.assembly
            }

        printt("making tracks and subtracks...")
        self._makeSubTracks()

    def _makeSubTracks(self):
        jobs = []
        for atn, labelNAndInfo in self.byAssayBiosampleType.iteritems():
            for labelN, info in labelNAndInfo.iteritems():
                jobs.append(merge_two_dicts(info,
                                            {"lookupByExp": self.lookupByExp,
                                             "idx": len(jobs) + 1,
                                             "total": len(self.byAssayBiosampleType)}))

        Parallel(n_jobs=self.args.j)(delayed(outputAllTracksByBiosampleType)
                                     (self.priority, job) for job in jobs)

    def _makeMainTrackDb(self):
        mainTrackDb = []

        for atn, labelNAndInfo in self.byAssayBiosampleType.iteritems():
            pri = self.priority.increment(1)
            totalExperiments = sum([len(info["exps"]) for info in labelNAndInfo.values()])
            shortLabel = self.labelNToNormal[atn]
            longLabel = self.labelNToNormal[atn] + " (%s experiments)" % totalExperiments

            mainTrackDb.append("""
track super_{atn}
superTrack on
priority {priority}
shortLabel {shortL}
longLabel {longL}
description {longL}
""".format(atn = atn,
           priority = pri,
           shortL=shortLabel,
           longL=Helpers.makeLongLabel(longLabel)))

        outF = StringIO.StringIO()
        outF.write('\n'.join(mainTrackDb))
        for atn, labelNAndInfo in self.byAssayBiosampleType.iteritems():
            for labelN, info in labelNAndInfo.iteritems():
                with open(info["fnp"]) as f:
                    outF.write(f.read())
                    outF.write('\n')
        return outF.getvalue()

def outputAllTracksByBiosampleType(priority, info):
    subGroups = outputSubTrack(priority, **info)
    info["subGroups"] = subGroups
    outputCompositeTrackByBiosampleType(**info)

def outputCompositeTrackByBiosampleType(assembly, assay_term_name, atn, label, labelN,
                                        exps, fnp, idx, total, subGroups, lookupByExp):
    if not os.path.exists(fnp):
        raise Exception("missing " + fnp)

    subGroupsDict = {}
    for k in Helpers.SubGroupKeys:
        subGroupsDict[k] = {a[0]:a[1] for a in subGroups[k]}
    longLabel = label + " (%s experiments)" % len(exps)

    print(subGroupsDict["assay"])

    subGroup1key = "biosample"
    subGroup2key = "age_sex"
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
    #     print("active biosample (composite):", labelNn)

    with open(fnp, 'w') as f:
        f.write("""
track {atn}_{labelN}
parent super_{atn}
compositeTrack on
centerLabelsDense on
""".format(atn=atn,
           labelN=labelN))
        if isActive:
            f.write("visibility full\n")
        f.write("""shortLabel {shortL}
description {longL}
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
""".format(shortL=Helpers.makeShortLabel(label),
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

def outputSubTrack(priority, assembly, assay_term_name, atn, label, labelN,
                   exps, fnp, idx, total, lookupByExp):
    actives = []
    # for expID in expIDs:
    #     if expID in lookupByExp:
    #         actives.append(lookupByExp[expID].isActive())
    isActive = any(t for t in actives)
    if isActive:
        print("active biosample:", btn)

    parent = Parent(atn + '_' + labelN, isActive)

    tracks = Tracks(assembly, parent, False)
    for exp in exps:
        active = False
        expID = exp.encodeID
        cREs = {}
        if expID in lookupByExp:
            active = lookupByExp[expID].isActive()
            cREs = lookupByExp[expID].cREs
        tracks.addExp(exp, True, cREs)

    Utils.ensureDir(fnp)
    with open(fnp, 'w') as f:
        for line in tracks.lines(priority):
            f.write(line)
    printWroteNumLines(fnp, idx, 'of', total)
    return tracks.subgroups()
