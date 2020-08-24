#!/usr/bin/env bash
computeOrderedMatrix -S ../scores.bw ../scores.bw -R ../regions.bed -o matrix_ordered.gz -g 1 --kmean 2
