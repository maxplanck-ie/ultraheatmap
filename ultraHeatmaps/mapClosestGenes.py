import pybedtools
from collections import defaultdict

def __splitClosestMapping(closestMapping, col_split):
    aInterval = [str(x).split('\t')[0:col_split] for x in closestMapping]

    bInterval = [str(x).split('\t')[col_split:] for x in closestMapping] #XXX what if there is no closest gene for a peak? the current code crashes if there is such a peak
    gff = [pybedtools.create_interval_from_list(i) for i in bInterval]

    return({'A': aInterval,'B': gff})

def __keymap_from_bed_and_gff(peaks, gff, key_definition = 'gene_id'):
    assert (len(peaks) == len(gff)), ("{} and {} are not the same size").format(len(peaks), len(gff))
    keyMap = defaultdict(lambda: [])
    for i in range(0,len(peaks)):
        ckey = ';'.join(peaks[i][0:7]) ##XXX shall we keep the number hard coded??
        cval = gff[i][key_definition]
        keyMap[ckey] = cval

    return(keyMap)

def keymap_from_closest_genes(closestMapping, peaks):
    assert (type(peaks) is type(pybedtools.BedTool())), ("{} is not class {}").format(type(peaks), type(pybedtools.BedTool()))
    splitDict = __splitClosestMapping(closestMapping, peaks.field_count())
    keyMap_closest = __keymap_from_bed_and_gff(splitDict['A'], splitDict['B'])
    return(keyMap_closest)
