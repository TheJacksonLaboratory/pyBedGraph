#!/bin/bash

sed -i -e 's/pyBigWig_exact /pyBW exact\n/g' $1/*
sed -i -e 's/pyBigWig_approx /pyBW app.\n/g' $1/*
sed -i -e 's/pyBigWig_approx/pyBW app./g' $1/*
sed -i -e 's/pyBedGraph_exact /pyBG exact\n/g' $1/*
sed -i -e 's/pBG approx./pyBG app./g' $1/*
