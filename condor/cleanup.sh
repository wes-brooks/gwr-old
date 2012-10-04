#! /bin/sh
for i in `ls output\\\*`; do
    mv $i "output/"`echo $i | awk -F \\\ {'print $2'}`;
done

