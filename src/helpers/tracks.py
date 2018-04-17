from __future__ import print_function

import sys
import os
import urllib
from collections import OrderedDict, defaultdict

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils'))
from utils import Utils, eprint, AddPath, printt, printWroteNumLines

import helpers as Helpers
from byAll import GetTissue, ColorByTissue

class LookupActive:
    def __init__(self, btid, btname, info, cREs):
        self.btid = btid
        self.btname = btname
        self.info = info
        self.cREs = cREs

    def isActive(self):
        return False

class LookupActiveForCcREs:
    def __init__(self, btid, btname, info, cREs):
        self.btid = btid
        self.btname = btname
        self.info = info
        self.cREs = cREs

    def isActive(self):
        r = self.btid in ["hepatocyte_derived_from_H9",
                          "bipolar_spindle_neuron_derived_from_induced_pluripotent_stem_cell",
                          "B_cell_adult"]
        return r

def outputLines(d, indentLevel, extras = {}):
    prefix = '\t' * indentLevel
    for k, v in d.iteritems():
        if v:
            yield prefix + k + " " + str(v) + '\n'
    for k, v in extras.iteritems():
        if v:
            yield prefix + k + " " + str(v) + '\n'
    yield '\n'

class Parent:
    def __init__(self, parent, on):
        self.parent = parent
        self.on = on

    def param(self, active):
        if active:
            return self.parent + ' on'
        return self.parent + ' off'

    def initials(self):
        return self.parent[:3] + '_'

class BigWigTrack(object):
    def __init__(self, assembly, exp, f, parent, active):
        self.assembly = assembly
        self.exp = exp
        self.f = f
        self.parent = parent
        self.active = active
        self.view = "bigWig"
        self.presentation = {}
        self.p = self._init()

    def _init(self):
        p = OrderedDict()
        p["track"] = self.parent.initials() + Helpers.sanitize(self.f.expID + '_' + self.f.fileID)
        p["parent"] = self.parent.param(self.active)
        p["subGroups"] = Helpers.unrollEquals(self._subgroups())
        p["bigDataUrl"] = self._url()
        p["visibility"] = Helpers.viz("full", self.active)
        p["type"] = "bigWig"
        p["color"] = Helpers.colorize(self.exp)
        p["height"] = "maxHeightPixels 64:12:8"
        p["shortLabel"] = Helpers.makeShortLabel(self.exp.assay_term_name, self.exp.tf)
        p["longLabel"] = Helpers.makeLongLabel(self.exp.assay_term_name + ' ' + self._desc())
        p["itemRgb"] = "On"
        p["darkerLabels"] = "on"
        p["metadata"] = Helpers.unrollEquals(self._metadata())
        p["view"] = self.view
        return p

    def _metadata(self):
        s = {}
        s["age"] = Helpers.sanitize(Helpers.getOrUnknown(self.exp.age_display))
        s["sex"] = Helpers.getOrUnknown(self.exp.donor_sex)
        s["accession"] = self.exp.encodeID
        s["description"] = Helpers.sanitize(self._desc())
        s["donor"] = self.exp.donor_id
        s["view"] = self.view
        return s

    def _subgroups(self):
        assay = self.exp.assay_term_name
        if "RNA-seq" == assay:
            assay = self.exp.assay_title
        target_label = ' '.join([self.exp.assay_term_name, self.exp.target, self.exp.label]).strip()
        s = {}
        s["donor"] = Helpers.getOrUnknown(self.exp.donor_id)
        s["assay"] = Helpers.getOrUnknown(assay)
        s["label"] = Helpers.getOrUnknown(self.exp.tf)
        s["target_label"] = Helpers.getOrUnknown(target_label)
        s["biosample"] = Helpers.getOrUnknown(self.exp.biosample_term_name)
        s["biosample_summary"] = Helpers.getOrUnknown(self.exp.biosample_summary).encode('ascii', 'ignore').decode('ascii')
        s["age"] = 'a' + Helpers.sanitize(Helpers.getOrUnknown(self.exp.age_display))
        s["sex"] = Helpers.getOrUnknown(self.exp.donor_sex)
        age_sex = ' '.join([e for e in [self.exp.age_display, self.exp.donor_sex] if e]).strip()
        s["age_sex"] = Helpers.getOrUnknown(age_sex)
        s["view"] = self.view
        self.presentation["label"] = (s["label"],
                                   Helpers.html_escape(Helpers.getOrUnknown(self.exp.tf)))
        self.presentation["assay"] = (s["assay"], s["assay"])
        self.presentation["donor"] = (s["donor"], s["donor"])
        self.presentation["target_label"] = (s["target_label"], s["target_label"])
        self.presentation["age"] = (s["age"],
                                    Helpers.html_escape(Helpers.getOrUnknown(self.exp.age_display)))
        self.presentation["view"] = (s["view"], s["view"])
        self.presentation["sex"] = (s["sex"], s["sex"])
        self.presentation["age_sex"] = (s["age_sex"], s["age_sex"])
        self.presentation["biosample"] = (s["biosample"], s["biosample"])
        self.presentation["biosample_summary"] = (s["biosample_summary"], s["biosample_summary"])
        self.presentation["tissue"] = self.presentation["biosample"]
        return s

    def _url(self):
        u = self.f.url
        if 'www.encodeproject.org' in u:
            if not u.endswith("?proxy=true"):
                u += "?proxy=true"
        return u

    def _desc(self):
        exp = self.exp
        desc = [self.exp.encodeID]
        if exp.biosample_summary:
            desc.append(Helpers.sanitize(exp.biosample_summary.strip()))
        elif exp.description:
            desc.append(exp.description)
        else:
            desc.append(exp.assay_term_name)
            if exp.tf:
                desc.append(exp.tf)
            age = exp.age_display
            if age and "unknown" != age:
                desc += [age]
        desc.append('(%s)' % self.f.output_type)
        return " ".join(desc)

    def lines(self, idx):
        extras = {}
        if self.active:
            extras["priority"] = idx
        return outputLines(self.p, 1, extras)

class BigWigTrackAll(BigWigTrack):
    def __init__(self, assembly, exp, f, parent, active, tissue):
        BigWigTrack.__init__(self, assembly, exp, f, parent, active)

        self.p["color"] = ColorByTissue(tissue)
        self.p["track"] = "all_" + self.p["track"]
        self.p["height"] = "maxHeightPixels 32:12:8"
        self.p["shortLabel"] = Helpers.makeShortLabel(tissue)

        self.presentation["tissue"] = (tissue, tissue)

class BigBedTrack(object):
    def __init__(self, assembly, exp, f, parent, active):
        self.assembly = assembly
        self.exp = exp
        self.f = f
        self.parent = parent
        self.active = active
        self.p = self._init()

    def _init(self):
        p = OrderedDict()
        p["track"] = self.parent.initials() + Helpers.sanitize(self.f.expID + '_' + self.f.fileID)
        p["parent"] = self.parent.param(self.parent.on)
        p["subGroups"] = Helpers.unrollEquals(self._subgroups())
        p["bigDataUrl"] = self._url()
        p["visibility"] = Helpers.viz("dense", self.active)
        p["type"] = "bigBed"
        p["shortLabel"] = Helpers.makeShortLabel(self.exp.assay_term_name, self.exp.tf)
        p["longLabel"] = Helpers.makeLongLabel(self._desc())
        p["itemRgb"] = "On"
        p["color"] = Helpers.colorize(self.exp)
        p["darkerLabels"] = "on"
        p["metadata"] = Helpers.unrollEquals(self._metadata())
        p["view"] = self.exp.encodeID
        return p

    def _url(self):
        u = self.f.url
        if 'www.encodeproject.org' in u:
            if not u.endswith("?proxy=true"):
                u += "?proxy=true"
        return u

    def _desc(self):
        exp = self.exp
        desc = [self.exp.encodeID]
        if exp.biosample_summary:
            desc.append(Helpers.sanitize(exp.biosample_summary.strip()))
        elif exp.description:
            desc.append(exp.description)
        else:
            desc.append(exp.assay_term_name)
            if exp.tf:
                desc.append(exp.tf)
            age = exp.age_display
            if age and "unknown" != age:
                desc += [age]
        desc.append('(%s)' % self.f.output_type)
        return " ".join(desc)

    def _metadata(self):
        s = {}
        s["age"] = Helpers.sanitize(Helpers.getOrUnknown(self.exp.age_display))
        s["sex"] = Helpers.getOrUnknown(self.exp.donor_sex)
        s["accession"] = self.exp.encodeID
        s["description"] = Helpers.sanitize(self._desc())
        s["donor"] = self.exp.donor_id
        s["view"] = self.exp.encodeID
        return s

    def _subgroups(self):
        assay = self.exp.assay_term_name
        if "RNA-seq" == assay:
            assay = self.exp.assay_title
        target_label = ' '.join([self.exp.assay_term_name, self.exp.target, self.exp.label]).strip()
        s = {}
        s["donor"] = Helpers.getOrUnknown(self.exp.donor_id)
        s["assay"] = Helpers.getOrUnknown(assay)
        s["label"] = Helpers.getOrUnknown(self.exp.tf)
        s["target_label"] = Helpers.getOrUnknown(target_label)
        s["biosample"] = Helpers.getOrUnknown(self.exp.biosample_term_name)
        s["biosample_summary"] = Helpers.getOrUnknown(self.exp.biosample_summary).encode('ascii', 'ignore').decode('ascii')
        s["age"] = 'a' + Helpers.sanitize(Helpers.getOrUnknown(self.exp.age_display))
        s["sex"] = Helpers.getOrUnknown(self.exp.donor_sex)
        age_sex = ' '.join([e for e in [self.exp.age_display, self.exp.donor_sex] if e]).strip()
        s["age_sex"] = Helpers.getOrUnknown(age_sex)
        s["view"] = self.exp.encodeID
        self.presentation = {}
        self.presentation["label"] = (s["label"],
                                   Helpers.html_escape(Helpers.getOrUnknown(self.exp.tf)))
        self.presentation["assay"] = (s["assay"], s["assay"])
        self.presentation["target_label"] = (s["target_label"], s["target_label"])
        self.presentation["donor"] = (s["donor"], s["donor"])
        self.presentation["age"] = (s["age"],
                                    Helpers.html_escape(Helpers.getOrUnknown(self.exp.age_display)))
        self.presentation["view"] = (s["view"], s["view"])
        self.presentation["sex"] = (s["sex"], s["sex"])
        self.presentation["age_sex"] = (s["age_sex"], s["age_sex"])
        self.presentation["biosample"] = (s["biosample"], s["biosample"])
        self.presentation["biosample_summary"] = (s["biosample_summary"], s["biosample_summary"])
        self.presentation["tissue"] = self.presentation["biosample"]
        return s

    def lines(self, idx):
        extras = {}
        if self.active:
            extras["priority"] = idx
        return outputLines(self.p, 2, extras)

class cRETrack(object):
    def __init__(self, assembly, exp, stateType, cREaccession, parent, active):
        self.assembly = assembly
        self.exp = exp
        self.stateType = stateType
        self.cREaccession = cREaccession
        self.parent = parent
        a = False
        if active and "5group" == stateType:
            a = True
        self.active = a
        self.p = self._init()

    def _init(self):
        p = OrderedDict()
        p["track"] = self.parent.initials() + Helpers.sanitize(self.exp.encodeID + '_' + self.cREaccession)
        p["parent"] = self.parent.param(self.parent.on)
        p["subGroups"] = Helpers.unrollEquals(self._subgroups())
        p["bigDataUrl"] = self._url()
        p["visibility"] = Helpers.viz("dense", self.active)
        p["type"] = "bigBed 9"
        p["shortLabel"] = Helpers.makeShortLabel(self.exp.assay_term_name, self.exp.tf)
        p["longLabel"] = Helpers.makeLongLabel(self._desc())
        p["itemRgb"] = "On"
        p["darkerLabels"] = "on"
        p["metadata"] = Helpers.unrollEquals(self._metadata())
        p["view"] = self.exp.encodeID
        return p

    def _url(self):
        return os.path.join("https://www.encodeproject.org/files/",
                            self.cREaccession,
                            "@@download/",
                            self.cREaccession + ".bigBed?proxy=true")

    def _desc(self):
        exp = self.exp
        return self.cREaccession + " " + self.stateType + " " + exp.description

    def _metadata(self):
        s = {}
        s["age"] = Helpers.sanitize(Helpers.getOrUnknown(self.exp.age_display))
        s["sex"] = Helpers.getOrUnknown(self.exp.donor_sex)
        s["accession"] = self.exp.encodeID
        s["description"] = Helpers.sanitize(self._desc())
        s["donor"] = self.exp.donor_id
        s["view"] = self.exp.encodeID
        return s

    def _subgroups(self):
        s = {}
        s["donor"] = Helpers.getOrUnknown(self.exp.donor_id)
        s["assay"] = Helpers.getOrUnknown(self.stateType)
        s["label"] = Helpers.getOrUnknown(self.exp.tf)
        s["biosample"] = Helpers.getOrUnknown(self.exp.biosample_term_name)
        s["age"] = 'a' + Helpers.sanitize(Helpers.getOrUnknown(self.exp.age_display))
        s["view"] = self.exp.encodeID
        self.presentation = {}
        self.presentation["label"] = (s["label"],
                                   Helpers.html_escape(Helpers.getOrUnknown(self.exp.tf)))
        self.presentation["assay"] = (s["assay"], s["assay"])
        self.presentation["donor"] = (s["donor"], s["donor"])
        self.presentation["age"] = (s["age"],
                                    Helpers.html_escape(Helpers.getOrUnknown(self.exp.age_display)))
        self.presentation["view"] = (s["view"], s["view"])
        self.presentation["biosample"] = (s["biosample"], s["biosample"])
        self.presentation["sex"] = ('', '')
        self.presentation["age_sex"] = ('', '')
        self.presentation["target_label"] = (s["assay"], s["assay"])
        self.presentation["biosample_summary"] = (s["biosample"], s["biosample"])
        self.presentation["tissue"] = self.presentation["biosample"]
        return s

    def lines(self, idx):
        extras = {}
        if self.active:
            extras["priority"] = idx
        return outputLines(self.p, 2, extras)

class CompositeExpTrack(object):
    def __init__(self, assembly, parent, exp, active):
        self.assembly = assembly
        self.parent = parent
        self.exp = exp
        self.active = active
        self.bedParent = Parent(parent.parent + '_view_' + exp.encodeID, parent.on)
        self.beds = []
        self.bigWigs = []
        self.cREs = []

    def _addExpBestBigWig(self, exp, active):
        files = Helpers.bigWigFilters(self.assembly, exp)
        expID = exp.encodeID

        ret = []
        if not files:
            eprint("missing bigwig for", expID)
        else:
            for f in files:
                t = BigWigTrack(self.assembly, exp, f, self.parent, active)
                ret.append(t)
        return ret

    def _addExpBestBigWigAll(self, exp, active):
        files = Helpers.bigWigFilters(self.assembly, exp)
        expID = exp.encodeID
        self.tissue = GetTissue(self.assembly, exp)

        ret = []
        if not files:
            eprint("missing bigwig for", expID)
        else:
            for f in files:
                t = BigWigTrackAll(self.assembly, exp, f, self.parent, active,
                                   self.tissue)
                ret.append(t)
        return ret

    def _addExpBestBed(self, exp, active):
        files = Helpers.bigBedFilters(self.assembly, exp)
        expID = exp.encodeID

        if not files:
            eprint("missing bed for", expID)
            #raise Exception("expected a file...", exp)
            return []
        ret = []
        for f in files:
            t = BigBedTrack(self.assembly, exp, f, self.bedParent, active)
            ret.append(t)
            break # TODO allow multiple bigBeds, if needed
        return ret

    def _addcREs(self, exp, active, cREs):
        ret = []
        cREaccessions = set()
        for stateType, accession in cREs.iteritems():
            if accession not in cREaccessions:
                t = cRETrack(self.assembly, exp, stateType, accession, self.bedParent, active)
                ret.append(t)
                cREaccessions.add(accession)
        return ret

    def addExp(self, cREs):
        self.beds = self._addExpBestBed(self.exp, self.active)
        self.bigWigs = self._addExpBestBigWig(self.exp, self.active)
        self.cREs = self._addcREs(self.exp, self.active, cREs)

    def addExpAll(self, cREs):
        self.bigWigs = self._addExpBestBigWigAll(self.exp, self.active)

    def view(self):
        p = OrderedDict()
        p["track"] = self.bedParent.parent
        p["parent"] = self.parent.param(self.active)
        p["view"] = self.exp.encodeID
        p["visibility"] = "dense"
        p["type"] = "bigBed"
        return p

    def tracks(self):
        for t in self.beds + self.cREs + self.bigWigs:
            yield t

class Tracks(object):
    def __init__(self, assembly, parent, priorityStart, isAll = False):
        self.assembly = assembly
        self.parent = parent
        self.priorityStart = priorityStart
        self.tracks = []
        self.isAll = isAll

    def addExp(self, exp, active, cREs):
        ct = CompositeExpTrack(self.assembly, self.parent, exp, active)
        ct.addExp(cREs)
        self.tracks.append(ct)

    def addExpAll(self, exp, active, cREs):
        ct = CompositeExpTrack(self.assembly, self.parent, exp, active)
        ct.addExpAll(cREs)
        self.tracks.append(ct)

    def lines(self):
        if self.isAll:
            tracks = self._sortAllTracks()
        else:
            tracks = self._sortedTracks()
        counter = 0
        for ct in tracks:
            for t in ct.bigWigs:
                counter += 1
                for line in t.lines(self.priorityStart + counter):
                    yield line
            if len(ct.beds + ct.cREs) > 0:
                # empty view not allowed
                for line in outputLines(ct.view(), 1):
                    yield line
                for t in ct.beds + ct.cREs:
                    counter += 1
                    for line in t.lines(self.priorityStart + counter):
                        yield line

    def _sortAllTracks(self):
        tracks = self.tracks

        def preferredSortOrder(track):
            return (track.tissue, track.exp.biosample_term_name)

        return sorted(tracks, key = lambda t: preferredSortOrder(t))

    def _sortedTracks(self):
        tracks = self.tracks

        def preferredSortOrder(exp):
            if exp.isDNaseSeq():
                return 1
            if exp.isChipSeqTF():
                if "CTCF" == exp.label:
                    return 4
            if exp.isChipSeqHistoneMark():
                if "H3K4me3" == exp.label:
                    return 2
                if "H3K27ac" == exp.label:
                    return 4
            return exp.label

        return sorted(tracks, key = lambda t: preferredSortOrder(t.exp))

    def subgroups(self):
        r = defaultdict(set)
        for tracks in self.tracks:
            for t in tracks.tracks():
                for k in Helpers.SubGroupKeys:
                    r[k].add(t.presentation[k])
        return r
