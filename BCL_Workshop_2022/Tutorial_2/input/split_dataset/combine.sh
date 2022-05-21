#!/bin/bash

# Don't assume they have zcat
gunzip 1843_*.sdf.gz
cat 1843_*.sdf > 1843_combined.sdf
gzip 1843_*.sdf
