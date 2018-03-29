import sys
import unittest

from gffannotator.GffAnnotator import *

class GffAnnotatorTest(unittest.TestCase):

    def test_class_call(self):
        genome='dm3'
        annotation='release78'
        gff_file = '../testing/data/dm3_ensembl78_chr4.gtf'
        FeatureDatabase = GffAnnotator(gff_file, genome, annotation, fast = True)
        dba_file = FeatureDatabase.db_file
        del FeatureDatabase
        self.assertTrue(not os.path.exists(dba_file))


    def test_class_functions(self):
        genome='dm3'
        annotation='release78'
        gff_file = '../testing/data/dm3_ensembl78_chr4.gtf'
        FeatureDatabase = GffAnnotator(gff_file, genome, annotation, fast = True)

        geneset = ['FBgn0040037','FBgn0052011','FBgn0052010','FBgn0017545']
        FeatureDatabase.geneId2Coordinates(geneset[:1])
        FeatureDatabase.geneId2Coordinates(geneset)

if __name__ == '__main__':
    unittest.main()
