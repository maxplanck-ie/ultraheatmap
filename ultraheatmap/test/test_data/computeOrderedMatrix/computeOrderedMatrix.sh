#!/usr/bin/env bash
computeOrderedMatrix -S ultraheatmap/test/test_data/scores.bw ultraheatmap/test/test_data/scores.bw -R ultraheatmap/test/test_data/regions.bed -o matrix_ordered.gz -g 1 --kmean 2
