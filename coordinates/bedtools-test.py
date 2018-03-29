import sys
import unittest

from shared.bedtools import *

class BedtoolsTest(unittest.TestCase):

    def test_parseGeneNames(self):
        tab1 = '../testing/data/DEseq2_dm3_chr4.tabular'
        gff = '../testing/data/dm3_ensembl78_chr4.gtf'
        deseq_tab = parseGeneIdTable(tab1)
        geneidList = deseq_tab['GeneID']
        bed = genes2Coordinates(geneidList[:5], gff, 'dm3', 'release78')
        for x in bed:
            print(x)


if __name__ == '__main__':
    unittest.main()
