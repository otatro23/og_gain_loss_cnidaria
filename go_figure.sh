#!/bin/bash

for dir in go_figure_tsvs*/* ; do 
	echo $dir
	gofigure.py -i $dir/BP_gain.tsv -j standard-plus -o $dir/gain -n bpo -e 100 -s user -su True
    gofigure.py -i $dir/BP_loss.tsv -j standard-plus -o $dir/loss -n bpo -e 100 -s user -su True
	gofigure.py -i $dir/CC_gain.tsv -j standard-plus -o $dir/gain -n cco -e 100 -s user -su True
	gofigure.py -i $dir/CC_loss.tsv -j standard-plus -o $dir/loss -n cco -e 100 -s user -su True
	gofigure.py -i $dir/MF_gain.tsv -j standard-plus -o $dir/gain -n mfo -e 100 -s user -su True
	gofigure.py -i $dir/MF_loss.tsv -j standard-plus -o $dir/loss -n mfo -e 100 -s user -su True
done

