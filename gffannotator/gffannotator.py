import gffutils
import os

class GffAnnotator:

    def __init__(self, gff_file, genome, annotation_version, fast = True):
        self.gff_file = gff_file
        self.db_file = ''.join([genome+'_', annotation_version, '.dba'])

        self.featureDb = None
        self.__createFeatureDatabase(fast = fast)

    ## should it really be treated as a temporary file?
    def __del__(self):
        os.remove(self.db_file)

    def __createFeatureDatabase(self, fast):
        if fast:
            self.__createDatabase(disable_infer_transcripts=True, disable_infer_genes=True)
        else:
            self.__createDatabase(disable_infer_transcripts=False, disable_infer_genes=False)

    def __createDatabase(self, disable_infer_transcripts, disable_infer_genes):
        try:
            db = gffutils.create_db(self.gff_file, dbfn=self.db_file,
                                keep_order=True, force = True,
                                disable_infer_transcripts=disable_infer_transcripts,
                                disable_infer_genes = disable_infer_genes)
        except: ## capture exceptions from gffutils properly a) wrong data b) already exists
            pass
        self.featureDb = gffutils.FeatureDB(dbfn=self.db_file, keep_order=True)

    def __geneid2Coord(self, geneid):
        return(self.featureDb[geneid])

    def geneId2Coordinates(self, geneids):
        return([self.__geneid2Coord(x) for x in geneids])

    def filter(self):
        pass

    def exportBedtools():
        pass
