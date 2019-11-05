#!/bin/bash
cd data
cat urls.txt | while read url ;
do
	wget -nc http://datasets.pacb.com.s3.amazonaws.com/2014/Iso-seq_Human_Tissues/${url}
done
mkdir raw 
mv m* raw
