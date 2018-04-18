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

from helpers.tracks import Tracks, Parent
import helpers.helpers as Helpers
from paths import Host, BaseWwwDir, BaseWwwTmpDir

sys.path.append(os.path.join(os.path.dirname(__file__), '../metadata/utils'))
from files_and_paths import Dirs
from utils import Utils, eprint, AddPath, printt, printWroteNumLines, dotdict
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

ActiveBiosamples = ["hepatocyte_derived_from_H9",
                    "bipolar_spindle_neuron_derived_from_induced_pluripotent_stem_cell",
                    "B_cell_adult"]

class MockFile:
    def __init__(self, eInfo, assembly):
        self.expID = eInfo["expID"]
        self.fileID = eInfo["fileID"]
        self.assembly = assembly
        self.output_type = "fold change over control"
        self.isPooled = True
        self.url = "https://www.encodeproject.org/files/{fileID}/@@download/{fileID}.bigWig".format(fileID = self.fileID)

    def isBigBed(self):
        return False

    def isBigWig(self):
        return True

    def isReleased(self):
        return True

class MockExp:
    def __init__(self, eInfo, assembly):
        for k, v in eInfo.items():
            setattr(self, k, v)
        self.assay_term_name = eInfo["assay"]
        self.encodeID = eInfo["expID"]
        self.files = [MockFile(eInfo, assembly)]
        self.donor_id = self.encodeID
        self.tf = self.assay
        self.label = self.assay
        self.target = self.assay
        self.age_display = ""
        self.donor_sex = ""
        self.description = eInfo["cellTypeDesc"]
        self.biosample_summary = eInfo["biosample_summary"]
        self.biosample_term_name = eInfo["cellTypeDesc"]
        self.ccREbigBeds = {}

    def isRnaSeqLike(self):
        return self.assay == "RNA-seq"

    def isDNaseSeq(self):
        return self.assay == "DNase-seq"

    def isChipSeq(self):
        return self.assay == "ChIP-seq"

    def isChipSeqTF(self):
        return self.assay == "CTCF"

    def isChipSeqHistoneMark(self):
        return self.assay == "H3K4me3" or self.assay == "H3K27ac"

def ccREexps(globalData, mw, assembly):
    creBigBeds = globalData["creBigBedsByCellType"]
    by4exps = globalData["byCellType"]

    expIDs = set()

    ret = []
    for ctn, eInfos in by4exps.iteritems():
        ctnExps = []
        for eInfo in eInfos:
            expID = eInfo["expID"]
            if expID in expIDs:
                eprint("skipping", expID)
                continue
            expIDs.add(expID) # b/c of ROADMAP
            e = MockExp(eInfo, assembly)
            e.active = ctn in ActiveBiosamples
            ctnExps.append(e)

        if not ctnExps:
            eprint("missing exps for", ctn)
            continue

        ccREbigBeds = creBigBeds.get(ctn, {})
        if not ccREbigBeds:
            eprint("missing ccREs for", ctn)
        ctnExps[0].ccREbigBeds = ccREbigBeds
        ret += ctnExps
    return ret

class TrackhubDbByCcREs:
    def __init__(self, args, assembly, globalData, mw):
        self.args = args
        self.assembly = assembly
        self.globalData = globalData
        self.mw = mw

        self.expsByAssay= [("candidate cis-Regulatory Regions", "ccres",
                            ccREexps),
        ]

        # assay x biosamepleType x biosamplesView

        self.byAssayBiosampleType = defaultdict(lambda: defaultdict(dict))
        self.subGroups = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.btToNormal = {}

    def run(self):
        for title, assayAbbr, expsF in self.expsByAssay:
            exps = expsF(self.globalData, self.mw, self.assembly)
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
                print(atn, bt)
                jobs.append(merge_two_dicts(info,
                                            {"idx": len(jobs) + 1,
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

def outputAllTracksByBiosampleType(info):
    subGroups = outputSubTrack(**info)
    info["subGroups"] = subGroups
    outputCompositeTrackByBiosampleType(**info)

def outputCompositeTrackByBiosampleType(assembly, assay_term_name,
                                        atn, biosample_type, bt,
                                        exps, fnp, idx, total,
                                        subGroups):
    if not os.path.exists(fnp):
        raise Exception("missing " + fnp)

    subGroupsDict = {}
    for k in Helpers.SubGroupKeys:
        subGroupsDict[k] = {a[0]:a[1] for a in subGroups[k]}
    longLabel = biosample_type + " (%s experiments)" % len(exps)

    if 0:
        for k, v in subGroupsDict.iteritems():
            print(k, v)

    subGroup1key = "biosample"
    subGroup2key = "assay"
    subGroup3key = "view"
    subGroup1 = Helpers.unrollEquals(subGroupsDict[subGroup1key])
    subGroup2 = Helpers.unrollEquals(subGroupsDict[subGroup2key])
    subGroup3 = Helpers.unrollEquals(subGroupsDict[subGroup3key])

    actives = []
    isActive = bt in ActiveBiosamples
    if isActive:
        print("active biosample (composite):", bt)

    with open(fnp) as f:
        subtracks = f.read()

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

def outputSubTrack(assembly, assay_term_name, atn, biosample_type,
                   bt, exps, fnp, idx, total):
    isActive = bt in ActiveBiosamples
    if isActive:
        print("active biosample:", bt)

    parent = Parent(atn + '_' + bt, isActive)

    tracks = Tracks(assembly, parent, (1 + idx) * 1000)
    for exp in exps:
        tracks.addExp(exp, exp.active, exp.ccREbigBeds)

    Utils.ensureDir(fnp)
    with open(fnp, 'w') as f:
        for line in tracks.lines():
            f.write(line)
    printWroteNumLines(fnp, idx, 'of', total)
    return tracks.subgroups()
