#! /bin/sh
for i in `ls $1/output-$2-*`; do
    tar -xf $i -C $1;
done

