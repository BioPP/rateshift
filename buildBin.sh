#! /bin/sh
arch=`uname -m`
version=0.9.0-1

strip rateshift/bpprateshift
tar cvzf bpprateshift-${arch}-bin-static-${version}.tar.gz rateshift/bpprateshift

