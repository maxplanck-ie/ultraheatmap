import pybedtools
from collections import defaultdict

def __splitClosestMapping(closestMapping, col_split):
    aInterval = [str(x).split('\t')[0:col_split] for x in closestMapping]
    bInterval = [str(x).split('\t')[col_split:] for x in closestMapping]
    gff = [pybedtools.create_interval_from_list(i) for i in bInterval]

    return({'A': aInterval,'B': gff})

def __keymap_from_bed_and_gff(peaks, gff):
    assert (len(peaks) == len(gff)), ("{} and {} are not the same size").format(len(peaks), len(gff))

    keyMap = defaultdict(lambda: [])
    for i in range(0,len(gff)):
        ckey = gff[i]['gene_id']
        cval = '_'.join(peaks[i][0:3])
        keyMap[ckey].append(cval)

    return(keyMap)

def keymap_from_closest_genes(closestMapping, peaks):
    ## should be replace by global BedTool-instance of peaksfile
    ncol1 = len(str(peaks[0]).strip().split('\t'))

    splitDict = __splitClosestMapping(closestMapping, ncol1)
    keyMap_closest = __keymap_from_bed_and_gff(splitDict['A'], splitDict['B'])
    return(keyMap_closest)
