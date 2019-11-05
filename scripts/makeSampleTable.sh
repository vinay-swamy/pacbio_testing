#!/bin/bash
# to get urls.txt, copy paste the downloads page  http://datasets.pacb.com.s3.amazonaws.com/2014/Iso-seq_Human_Tissues/list.html

tmp=/tmp/asdfasfd.txt
grep raw urls.txt | grep 'bax' - | tr '/' '\t' | awk '{OFS="\t"};{print $2,$3}' - > $tmp
names=/tmp/asdfasdfasdfafd.txt
cut -f2 $tmp | tr '_' '\t'| awk '{OFS="_"};{print $1,$2}' - >$names
tissue=/tmp/dfgdsf.txt
cut -f1 $tmp > $tissue

path=/tmp/dsfgdfgsdfgsdfg.txt
cut -f2 $tmp > $path

paste $names $tissue $path 

