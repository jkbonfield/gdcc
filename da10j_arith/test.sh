#!/bin/bash

mkdir  ../T2_demo.enc ../T2_demo.dec 2>/dev/null
rm ../T2_demo.enc/* ../T2_demo.dec/* 2>/dev/null

time for i in ../T2_demo/*
do
    f=`echo $i | sed 's#.*T2_demo/##'`
    printf "%-40s " $i
    ./da9j e $i ../T2_demo.enc/$f 2>/dev/null
    cat ../T2_demo.enc/$f| wc -c
done

cat ../T2_demo.enc/*|wc -c

time for i in ../T2_demo.enc/*
do f=`echo $i | sed 's#.*T2_demo.enc/##'`
   ./da9j d $i ../T2_demo.dec/$f 2>/dev/null
done

cd ../T2_demo
for i in *
do
    cmp $i ../T2_demo.dec/$i
done
