## proovread Makefile
## Author: thackl@lim4.de

all: bwa string

bwa:
	cd util/bwa; make;

string:
	cd src/String-Similarity-1.04; \
	perl Makefile.PL INSTALLSITEARCH=../../lib INSTALLSITEMAN3DIR=../../lib/man; \
	make; make install;

clean:
	cd src/String-Similarity-1.04; \
	make clean;


