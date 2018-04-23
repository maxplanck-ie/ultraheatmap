import gffutils
from pybedtools import BedTool
from collections import defaultdict

import os
import sys

class GffAnnotator:
    def __init__(self, gff_file, fast = True, verbose = False):
        self.gff_file = gff_file
        print(self.gff_file)
        self.featureDb = None
        self.__createFeatureDatabase(fast = fast)

        self.filters = {}
        self.filters["featuretypes"] = [x for x in self.featureDb.featuretypes()]
        self.filters["gene_biotypes"] = set([x.attributes['gene_biotype'][0] for x in self.featureDb.all_features()])

    def __createFeatureDatabase(self, fast):
        if fast:
            self.__createDatabase(disable_infer_transcripts=True, disable_infer_genes=True)
        else:
            self.__createDatabase(disable_infer_transcripts=False, disable_infer_genes=False)

    def __createDatabase(self, disable_infer_transcripts, disable_infer_genes):
        try:
            self.featureDb = gffutils.create_db(self.gff_file, dbfn=":memory:",
                                keep_order=True, force = False,
                                disable_infer_transcripts=disable_infer_transcripts,
                                disable_infer_genes = disable_infer_genes)
        except:
            pass

    def __geneid2Coord(self, geneid):
        try:
            return(self.featureDb[geneid])
        except:
            print ("Warning: %s not found" %geneid, file=sys.stderr)

    def __bed12(self, feature, stream):
        try:
            return(self.featureDb.bed12(feature))
        except:
            print ("Warning: %s not found" % feature, file=sys.stderr)

    def geneId2Coordinates(self, geneids):
        return([self.__geneid2Coord(x) for x in geneids])

    def geneid2BedTool(self, geneIds, filename = None, as_pybedtool = False):
        featureSet = self.geneId2Coordinates(geneIds)
        featureCoords = BedTool('\n'.join([self.__bed12(feature, sys.stdout) for feature in featureSet]), from_string = True)
        if filename:
            sys.stderr.write('Writing to file:\n' + filename)
            with open(filename, 'w') as bed12:
                bed12.write(str(featureCoords))
        return(featureCoords)

    ## keymap keys:
    ## gene keys: gff gene_id
    ## peak keys: <chr>_<start>_<end>
    def geneid2keymap(self, geneids):
        keyMap = defaultdict(lambda: None)
        for gid in geneids:
            feature= self.__geneid2Coord(gid)
            keyMap[gid] = "{}_{}_{}".format(feature.seqid, feature.start, feature.end)
        return(keyMap)

    def filter(self):
        pass
