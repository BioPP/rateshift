#! /bin/sh
arch=`uname -m`
version=1.0.1-1

strip rateshift/bpprateshift
tar cvzf bpprateshift-${arch}-bin-static-${version}.tar.gz rateshift/bpprateshift

