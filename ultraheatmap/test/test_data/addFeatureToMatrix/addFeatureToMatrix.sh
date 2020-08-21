#!/usr/bin/env bash
computeMatrix reference-point -S scores.bw -R regions.bed -o matrix.gz
addFeatureToMatrix -m matrix.gz -o appended_matrix.gz -t features.bed -f score1 score2 --featureIdColumn name
