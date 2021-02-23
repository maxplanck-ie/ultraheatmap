import pybedtools
from collections import defaultdict

def __splitClosestMapping(closestMapping, col_split):
    aInterval = []
    bInterval = []

    for x in closestMapping:
        a = str(x).split('\t')[0:col_split]
        b = str(x).split('\t')[col_split:]
        if '-1' not in b:
            aInterval.append(';'.join(str(v) for v in a))
            bInterval.append(b)
        else:
            print("region ", a , "got no gene to be found as its closest gene.")

    gff = [pybedtools.create_interval_from_list(i) for i in bInterval]

    return({'A': aInterval,'B': gff})


def __keymap_from_bed_and_gff(peaks, aInterval, gff, key_definition = 'gene_id'):
    count = 0
    keyMap = defaultdict(lambda: [])

    for peak in peaks:
        ckey = ';'.join(str(v) for v in peak)
        if ckey in aInterval:
            i = aInterval.index(ckey)
            cval = gff[i][key_definition]
        else:
            cval = "no_gene"
            count += 1
        keyMap[ckey] = cval

    return(keyMap)


def keymap_from_closest_genes(closestMapping, peaks, field_count):
    splitDict = __splitClosestMapping(closestMapping, field_count)
    keyMap_closest = __keymap_from_bed_and_gff(peaks, splitDict['A'], splitDict['B'])
    return(keyMap_closest)
