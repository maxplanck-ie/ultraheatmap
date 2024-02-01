#!/usr/bin/env bash
computeMatrix reference-point -S ultraheatmap/test/test_data/scores.bw -R ultraheatmap/test/test_data/regions.bed -o matrix.gz
addFeatureToMatrix -m matrix.gz -o appended_matrix.gz -t ultraheatmap/test/test_data/addFeatureToMatrix/features.bed -f score1 score2 --featureIdColumn name
