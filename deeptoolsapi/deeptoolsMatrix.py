import numpy as np
import gzip
import sys
## full class copied from deeptools
# https://github.com/deeptools/deepTools/blob/master/deeptools/heatmapper.py
def toBytes(s):
        """
        Like toString, but for functions requiring bytes in python3
        """
        if sys.version_info[0] == 2:
            return s
        if isinstance(s, bytes):
            return s
        if isinstance(s, str):
            return bytes(s, 'ascii')
        if isinstance(s, list):
            return [toBytes(x) for x in s]
        return s

def toString(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')
    if isinstance(s, list):
        return [toString(x) for x in s]
    return s

def read_matrix_file(matrix_file):
    # reads a bed file containing the position
    # of genomic intervals
    # In case a hash sign '#' is found in the
    # file, this is considered as a delimiter
    # to split the heatmap into groups
    parameters = dict()
    special_params = set(['unscaled 5 prime', 'unscaled 3 prime', 'body', 'downstream', 'upstream', 'ref point', 'bin size'])
    import json
    regions = []
    matrix_rows = []
    current_group_index = 0
    max_group_bound = None
    fh = gzip.open(matrix_file)

    for line in fh:
        line = toString(line).strip()
        # read the header file containing the parameters
        # used
        if line.startswith("@"):
            # the parameters used are saved using
            # json
            parameters = json.loads(line[1:].strip())
            max_group_bound = parameters['group_boundaries'][1]
            continue
            # split the line into bed interval and matrix values
        region = line.split('\t')
        chrom, start, end, name, score, strand = region[0:6]
        matrix_row = np.ma.masked_invalid(np.fromiter(region[6:], np.float))
        matrix_rows.append(matrix_row)
        starts = start.split(",")
        ends = end.split(",")
        regs = [(int(x), int(y)) for x, y in zip(starts, ends)]
        # get the group index
        if len(regions) >= max_group_bound:
            current_group_index += 1
            max_group_bound = parameters['group_boundaries'][current_group_index + 1]
        regions.append([chrom, regs, name, max_group_bound, strand, score])
        matrix = np.vstack(matrix_rows)
    matrix = Matrix(regions, matrix, parameters['group_boundaries'],
                          parameters['sample_boundaries'],
                          group_labels=parameters['group_labels'],
                          sample_labels=parameters['sample_labels'])
                # Versions of computeMatrix before 3.0 didn't have an entry of these per column, fix that
    nSamples = len(matrix.sample_labels)
    h = dict()
    for k, v in parameters.items():
        if k in special_params and type(v) is not list:
            v = [v] * nSamples
            if len(v) == 0:
                v = [None] * nSamples
        h[k] = v
    print("h: "+ str(h))
    print(len(matrix.regions))
    return matrix, h


class Matrix:

    def __init__(self, regions, matrix, group_boundaries, sample_boundaries,
                 group_labels=None, sample_labels=None):
        # simple checks
        assert matrix.shape[0] == group_boundaries[-1], \
            "row max do not match matrix shape"
        assert matrix.shape[1] == sample_boundaries[-1], \
            "col max do not match matrix shape"
        assert matrix.shape[0] == len(regions)

        self.regions = regions
        self.matrix = matrix
        self.group_boundaries = group_boundaries # index vector
        self.sample_boundaries = sample_boundaries # index vector
        self.sort_method = None
        self.sort_using = None

        if group_labels is None:
            self.group_labels = ['group {}'.format(x)
                                 for x in range(len(group_boundaries) - 1)]
        else:
            assert len(group_labels) == len(group_boundaries) - 1, \
                "number of group labels does not match number of groups"
            self.group_labels = group_labels

        if sample_labels is None:
            self.sample_labels = ['sample {}'.format(x)
                                  for x in range(len(sample_boundaries) - 1)]
        else:
            assert len(sample_labels) == len(sample_boundaries) - 1, \
                "number of sample labels does not match number of samples"
            self.sample_labels = sample_labels


    def get_matrix(self, group, sample):
        """
        Returns a sub matrix from the large
        matrix. Group and sample are ids,
        thus, row = 0, col=0 get the first group
        of the first sample.
        Returns
        -------
        dictionary containing the matrix,
        the group label and the sample label
        """
        group_start = self.group_boundaries[group]
        group_end = self.group_boundaries[group + 1]
        sample_start = self.sample_boundaries[sample]
        sample_end = self.sample_boundaries[sample + 1]

        return {'matrix': np.ma.masked_invalid(self.matrix[group_start:group_end, :][:, sample_start:sample_end]),
                'group': self.group_labels[group],
                'sample': self.sample_labels[sample]}

    def get_num_samples(self):
        return len(self.sample_labels)

    def get_num_groups(self):
        return len(self.group_labels)

    def set_group_labels(self, new_labels):
        """ sets new labels for groups
        """
        if len(new_labels) != len(self.group_labels):
            raise ValueError("length new labels != length original labels")
        self.group_labels = new_labels

    def set_sample_labels(self, new_labels):
        """ sets new labels for groups
        """
        if len(new_labels) != len(self.sample_labels):
            raise ValueError("length new labels != length original labels")
        self.sample_labels = new_labels

    def set_sorting_method(self, sort_method, sort_using):
        self.sort_method = sort_method
        self.sort_using = sort_using

    def get_regions(self):
        """Returns the regions per group
        Returns
        ------
        list
            Each element of the list is itself a list
            of dictionaries containing the regions info:
            chrom, start, end, strand, name etc.
            Each element of the list corresponds to each
            of the groups
        """
        regions = []
        for idx in range(len(self.group_labels)):
            start = self.group_boundaries[idx]
            end = self.group_boundaries[idx + 1]
            regions.append(self.regions[start:end])

        return regions

    def sort_groups(self, sort_using='mean', sort_method='no', sample_list=None):
        """
        Sorts and rearranges the submatrices according to the
        sorting method given.
        """
        if sort_method == 'no':
            return

        if (sample_list is not None) and (len(sample_list) > 0):
            # get the ids that correspond to the selected sample list
            idx_to_keep = []
            for sample_idx in sample_list:
                idx_to_keep += range(self.sample_boundaries[sample_idx], self.sample_boundaries[sample_idx + 1])

            matrix = self.matrix[:, idx_to_keep]

        else:
            matrix = self.matrix

        # compute the row average:
        if sort_using == 'region_length':
            matrix_avgs = list()
            for x in self.regions:
                matrix_avgs.append(np.sum([bar[1] - bar[0] for bar in x[1]]))
            matrix_avgs = np.array(matrix_avgs)
        elif sort_using == 'mean':
            matrix_avgs = np.nanmean(matrix, axis=1)
        elif sort_using == 'mean':
            matrix_avgs = np.nanmean(matrix, axis=1)
        elif sort_using == 'median':
            matrix_avgs = np.nanmedian(matrix, axis=1)
        elif sort_using == 'max':
            matrix_avgs = np.nanmax(matrix, axis=1)
        elif sort_using == 'min':
            matrix_avgs = np.nanmin(matrix, axis=1)
        elif sort_using == 'sum':
            matrix_avgs = np.nansum(matrix, axis=1)
        else:
            sys.exit("{} is an unsupported sorting method".format(sort_using))

        # order per group
        _sorted_regions = []
        _sorted_matrix = []
        for idx in range(len(self.group_labels)):
            start = self.group_boundaries[idx]
            end = self.group_boundaries[idx + 1]
            order = matrix_avgs[start:end].argsort()
            if sort_method == 'descend':
                order = order[::-1]
            _sorted_matrix.append(self.matrix[start:end, :][order, :])
            # sort the regions
            _reg = self.regions[start:end]
            for idx in order:
                _sorted_regions.append(_reg[idx])

        self.matrix = np.vstack(_sorted_matrix)
        self.regions = _sorted_regions
        self.set_sorting_method(sort_method, sort_using)

    def hmcluster(self, k, method='kmeans'):

        matrix = np.asarray(self.matrix)
        if np.any(np.isnan(matrix)):
            # replace nans for 0 otherwise kmeans produces a weird behaviour
            sys.stderr.write("*Warning* For clustering nan values have to be replaced by zeros \n")
            matrix[np.isnan(matrix)] = 0

        if method == 'kmeans':
            from scipy.cluster.vq import vq, kmeans

            centroids, _ = kmeans(matrix, k)
            # order the centroids in an attempt to
            # get the same cluster order
            cluster_labels, _ = vq(matrix, centroids)

        if method == 'hierarchical':
            # normally too slow for large data sets
            from scipy.cluster.hierarchy import fcluster, linkage
            Z = linkage(matrix, method='ward', metric='euclidean')
            cluster_labels = fcluster(Z, k, criterion='maxclust')
            # hierarchical clustering labels from 1 .. k
            # while k-means labels 0 .. k -1
            # Thus, for consistency, we subtract 1
            cluster_labels -= 1

        # sort clusters
        _clustered_mean = []
        _cluster_ids_list = []
        for cluster in range(k):
            cluster_ids = np.flatnonzero(cluster_labels == cluster)
            _cluster_ids_list.append(cluster_ids)
            _clustered_mean.append(self.matrix[cluster_ids, :].mean())

        # reorder clusters based on mean
        cluster_order = np.argsort(_clustered_mean)[::-1]

        # create groups using the clustering
        self.group_labels = []
        self.group_boundaries = [0]
        _clustered_regions = []
        _clustered_matrix = []
        cluster_number = 1
        for cluster in cluster_order:
            self.group_labels.append("cluster_{}".format(cluster_number))
            cluster_number += 1
            cluster_ids = _cluster_ids_list[cluster]
            self.group_boundaries.append(self.group_boundaries[-1] +
                                         len(cluster_ids))
            _clustered_matrix.append(self.matrix[cluster_ids, :])
            for idx in cluster_ids:
                _clustered_regions.append(self.regions[idx])

        self.regions = _clustered_regions
        self.matrix = np.vstack(_clustered_matrix)
        return idx

    def removeempty(self):
        """
        removes matrix rows containing only zeros or nans
        """
        to_keep = []
        score_list = np.ma.masked_invalid(np.mean(self.matrix, axis=1))
        for idx, region in enumerate(self.regions):
            if np.ma.is_masked(score_list[idx]) or np.float(score_list[idx]) == 0:
                continue
            else:
                to_keep.append(idx)
        self.regions = [self.regions[x] for x in to_keep]
        self.matrix = self.matrix[to_keep, :]
        # adjust sample boundaries
        to_keep = np.array(to_keep)
        self.group_boundaries = [len(to_keep[to_keep < x]) for x in self.group_boundaries]

    def flatten(self):
        """
        flatten and remove nans from matrix. Useful
        to get max and mins from matrix.
        :return flattened matrix
        """
        matrix_flatten = np.asarray(self.matrix.flatten())
        # nans are removed from the flattened array
        matrix_flatten = matrix_flatten[~np.isnan(matrix_flatten)]
        if len(matrix_flatten) == 0:
            num_nan = len(np.flatnonzero(np.isnan(self.matrix.flatten())))
            raise ValueError("matrix only contains nans "
"(total nans: {})".format(num_nan))



    def save_matrix(self, file_name, mode = 'scale-regions', ref_point = None):
        """
        saves the data required to reconstruct the matrix
        the format is:
        A header containing the parameters used to create the matrix
        encoded as:
        @key:value\tkey2:value2 etc...
        The rest of the file has the same first 5 columns of a
        BED file: chromosome name, start, end, name, score and strand,
        all separated by tabs. After the fifth column the matrix
        values are appended separated by tabs.
        Groups are separated by adding a line starting with a hash (#)
        and followed by the group name.
        The file is gzipped.
        """
        import json
        parameters=dict()
        parameters['sample_labels'] = self.sample_labels
        parameters['group_labels'] = self.group_labels
        parameters['sample_boundaries'] = self.sample_boundaries
        parameters['group_boundaries'] = self.group_boundaries
#        special_params = set(['unscaled 5 prime', 'unscaled 3 prime', 'body', 'downstream', 'upstream', 'ref point', 'bin size'])
        special_params = {'unscaled 5 prime':0, 'unscaled 3 prime':0, 'body':1000, 'downstream':0, 'upstream':0, 'ref point':None, 'bin size':10}
        if mode == 'reference-point':
            special_params['ref point'] = ref_point
        nSamples = len(self.sample_labels)
        h = dict()
        for k, v in parameters.items():
            if type(v) is list and len(v) == 0:
                v = None
            h[k] = v
        for k , v in special_params.items():
            v = [v] * nSamples
            h[k] = v

        print(h)
        fh = gzip.open(file_name, 'wb')
        params_str = json.dumps(h, separators=(',', ':'))
        fh.write(toBytes("@" + params_str + "\n"))
        score_list = np.ma.masked_invalid(np.mean(self.matrix, axis=1))
        for idx, region in enumerate(self.regions):
            # join np_array values
            # keeping nans while converting them to strings
            if not np.ma.is_masked(score_list[idx]):
                np.float(score_list[idx])
            matrix_values = "\t".join(np.char.mod('%f', self.matrix[idx, :]))
            starts = ["{0}".format(x[0]) for x in region[1]]
            ends = ["{0}".format(x[1]) for x in region[1]]
            starts = ",".join(starts)
            ends = ",".join(ends)
            # BEDish format (we don't currently store the score)
            fh.write(
                toBytes('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                        region[0],
                        starts,
                        ends,
                        region[2],
                        region[5],
                        region[4],
                        matrix_values)))
        fh.close()
