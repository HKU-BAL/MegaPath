#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
	echo "Usage $0 <output_prefix>";
	exit 1;
fi

cat $1.dpout.1
egrep -v '^@' $1.unpair