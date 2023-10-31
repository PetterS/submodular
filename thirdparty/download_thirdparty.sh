#!/bin/bash

#wget http://pub.ist.ac.at/~vnk/software/maxflow-v3.01.src.tar.gz 
#tar -xzf maxflow-v3.01.src.tar.gz 
wget http://vision.csd.uwo.ca/files/maxflow-v3.01.zip
unzip maxflow-v3.01.zip -d maxflow-tmp
mv maxflow-tmp/* maxflow-v3.01.src/
rm -rf maxflow-tmp

wget http://pub.ist.ac.at/~vnk/software/QPBO-v1.3.src.tar.gz 
tar -xzf QPBO-v1.3.src.tar.gz 
mv QPBO-v1.3.src/* QPBO  
cp QPBO/instances.inc.new QPBO/instances.inc

wget http://www.f.waseda.jp/hfs/HOCR1.02.zip  
unzip HOCR1.02.zip -d tmp 
mv tmp/Image.h HOCR/Image.h 
mv tmp/HOCR/HOCR.h HOCR/HOCR.h 
mv tmp/HOCR/HOCR0.h HOCR/HOCR0.h 

rm QPBO-v1.3.src.tar.gz 
rm HOCR1.02.zip 
rm -r tmp 
rm -r QPBO-v1.3.src 

