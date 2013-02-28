#! /bin/sh
for i in 410 411 412 415 419 420 421 422 1638
do
tar -zvcf output-43-$i.tgz2 `find ./ -regex "^[^o].*\.43\.$i.*"`
done




