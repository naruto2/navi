#!/bin/sh

function testnavi() {
    echo -n $1
    diff $1.gz ../src/$1.gz > ../tmp.log
    if [ -s ../tmp.log ];then
	echo  "  ..... NG"
    else
	echo  "  ..... OK"
    fi
    rm -f  ../tmp.log
}


testnavi 00000
testnavi 00001
testnavi 00002
testnavi 00003
testnavi 00004
testnavi 00005
