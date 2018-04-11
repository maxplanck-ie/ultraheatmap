import gffutils
import os

class GffAnnotator:

    def __init__(self, gff_file, fast = True, verbose = False):
        self.gff_file = gff_file

        self.featureDb = None
        self.__createFeatureDatabase(fast = fast)
        self.featureTypes = [x for x in self.featureDb.featuretypes()]
        if verbose:
            print(self.featureTypes)

        self.filters = list()
        self.filters.featuretype = [x for x in self.featureDb.featureTypes()]
        self.filters.gene_biotypes = set([x.attributes['gene_biotype'][0] for x in self.featureDb.all_features()])


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
        return(self.featureDb[geneid])

    def geneId2Coordinates(self, geneids):
        return([self.__geneid2Coord(x) for x in geneids])

    def exportBed12(self, filename, geneids):
        with open(filename, 'w') as bed12:
            for gid in geneids:
                try:
                    bed12.write(self.featureDb.bed12(gid) + os.linesep)
                except:
                     print ("Warning: %s not found" % gid, file=sys.stderr)

    def filter(self):
        pass
