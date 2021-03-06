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

class TrackhubDbByAssayByBiosampleType:
    def __init__(self, args, assembly, globalData, mw, priority):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData
        self.mw = mw
        self.priority = priority

        self.expsByAssay= [("DNase-seq", "dnase",
                            "DNase-seq", True,
                            self.mw.dnases_useful),
                           ("Histone by Biosample", "histone_modifications",
                            "Histone modifications and variants", False,
                            self.mw.chipseq_histones_useful),
                           ("RNA-seq", "transcription",
                            "RNA-seq", True,
                            self.mw.transcription_useful),
                           ("microRNA-seq", "microRNAseq",
                            "microRNA-seq", True,
                            self.mw.microRNAseq_useful),
                           ("TFs by Biosample Type", "transcription_factors",
                            "Transcription Factors", False,
                            self.mw.chipseq_tfs_useful)
        ]
        if "mm10" == assembly:
            self.expsByAssay.append(("ATAC-seq", "atac_seq",
                                     "ATAC-seq", True,
                                     self.mw.atac_seq_useful))
        if "hg19" == assembly:
            self.expsByAssay.append(("RAMPAGE", "rampage",
                                     "RAMPAGE", True,
                                     self.mw.rampage_useful))


        # assay x biosamepleType x biosamplesView

        self.byAssayBiosampleType = defaultdict(lambda: defaultdict(dict))
        self.subGroups = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.btToNormal = {}
        self.lookupByExp = {}

    def run(self):
        for title, assayAbbr, longLabelBase, showAllTrack, expsF in self.expsByAssay:
            exps = expsF()
            self._build(title, assayAbbr, showAllTrack, exps, longLabelBase)
        return self._makeMainTrackDb()

    def _build(self, assay_term_name, atn, showAllTrack, exps, longLabelBase):
        printt("building", assay_term_name, "...")

        def sorter(exp):
            return (exp.biosample_type)
        exps.sort(key = sorter)

        self.btToNormal[atn] = assay_term_name

        if showAllTrack:
            biosample_type = "_ALL DATA"
            bt = "0_all"
            self.btToNormal[bt] = biosample_type

            fnp = os.path.join(BaseWwwTmpDir, self.assembly, "subtracks", atn, bt +'.txt')
            self.byAssayBiosampleType[atn][bt] = {
                "assay_term_name": assay_term_name,
                "atn": atn,
                "biosample_type": biosample_type,
                "bt": bt,
                "fnp": fnp,
                "exps": exps,
                "assembly": self.assembly,
                "longLabelBase": longLabelBase
            }

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
                "assembly": self.assembly,
                "longLabelBase": longLabelBase
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

        Parallel(n_jobs=self.args.j)(delayed(outputAllTracksByBiosampleType)
                                     (self.priority, job) for job in jobs)

    def _makeMainTrackDb(self):
        mainTrackDb = []

        for atn, btAndInfo in self.byAssayBiosampleType.iteritems():
            pri = self.priority.increment(1)

            totalExperiments = 0
            for bt, info in btAndInfo.iteritems():
                if "0_all" != bt:
                    totalExperiments += len(info["exps"])

            shortLabel = self.btToNormal[atn]
            longLabel = info["longLabelBase"] + " (%s experiments)" % totalExperiments

            mainTrackDb.append("""
track super_{atn}
superTrack on
priority {priority}
shortLabel {shortL}
longLabel {longL}
""".format(atn = atn,
           priority = pri,
           shortL=shortLabel,
           longL=Helpers.makeLongLabel(longLabel)))

        outF = StringIO.StringIO()
        outF.write('\n'.join(mainTrackDb))
        for atn, btAndInfo in self.byAssayBiosampleType.iteritems():
            for bt, info in btAndInfo.iteritems():
                with open(info["fnp"]) as f:
                    outF.write(f.read())
                    outF.write('\n')
        return outF.getvalue()

def outputAllTracksByBiosampleType(priority, info):
    subGroups = outputSubTrack(priority, **info)
    info["subGroups"] = subGroups
    outputCompositeTrackByBiosampleType(**info)

def outputCompositeTrackByBiosampleType(assembly, assay_term_name, atn, biosample_type, bt,
                                        exps, fnp, idx, total, subGroups, lookupByExp, longLabelBase = None):
    if not os.path.exists(fnp):
        raise Exception("missing " + fnp)

    subGroupsDict = {}
    for k in Helpers.SubGroupKeys:
        subGroupsDict[k] = {a[0]:a[1] for a in subGroups[k]}
    longLabel = biosample_type + " (%s experiments)" % len(exps)

    print(subGroupsDict["assay"])

    if atn in ["dnase", "atac_seq"]:
        subGroup1key = "biosample"
        subGroup2key = "age"
    elif "histone_modifications" == atn:
        subGroup1key = "biosample"
        subGroup2key = "label"
    elif "transcription_factors" == atn:
        subGroup1key = "biosample"
        subGroup2key = "label"
    elif atn in ["transcription", "rampage", "microRNAseq"]:
        subGroup1key = "biosample"
        subGroup2key = "assay"
    else:
        subGroup1key = "donor"
        subGroup2key = "age"

    if "0_all" == bt:
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
    #     print("active biosample (composite):", btn)

    with open(fnp, 'w') as f:
        f.write("""
track {atn}_{bt}
parent super_{atn}
compositeTrack on
centerLabelsDense on
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

def outputSubTrack(priority, assembly, assay_term_name, atn, biosample_type, bt,
                   exps, fnp, idx, total, lookupByExp, longLabelBase = None):
    actives = []
    # for expID in expIDs:
    #     if expID in lookupByExp:
    #         actives.append(lookupByExp[expID].isActive())
    isActive = any(t for t in actives)
    if isActive:
        print("active biosample:", btn)

    parent = Parent(atn + '_' + bt, isActive)

    tracks = Tracks(assembly, parent, "0_all" == bt)
    for exp in exps:
        active = False
        expID = exp.encodeID
        cREs = {}
        if expID in lookupByExp:
            active = lookupByExp[expID].isActive()
            cREs = lookupByExp[expID].cREs
        if "0_all" == bt:
            tracks.addExpAll(exp, True, cREs)
        else:
            tracks.addExp(exp, True, cREs)

    Utils.ensureDir(fnp)
    with open(fnp, 'w') as f:
        for line in tracks.lines(priority):
            f.write(line)
    printWroteNumLines(fnp, idx, 'of', total)
    return tracks.subgroups()
