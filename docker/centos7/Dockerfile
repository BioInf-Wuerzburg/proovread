FROM centos:centos7

# getting the last blast+ package
# install required packages
RUN yum install --assumeyes wget
RUN wget -O /tmp/blast.rpm ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$(wget http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST -O -| grep -vP 'src|md5' | grep -oP '(?<=<a href=")ncbi.*rpm') && yum install --assumeyes /tmp/blast.rpm

# getting the samtools 1.2
# install required packages
RUN yum install --assumeyes tar bzip2 ncurses-devel gcc make zlib-devel
RUN wget -O /tmp/samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
RUN cd /tmp/ && tar xjf samtools.tar.bz2 && cd samtools-1.2 && make && make install

# getting the last proovread releas
RUN yum install --assumeyes git cpan which perl-Time-HiRes perl-File-Which
RUN perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit' && cpan Log::Log4perl
RUN cd ~ && git clone --recursive https://github.com/BioInf-Wuerzburg/proovread && cd proovread && git checkout proovread-2.12 && cd util/bwa && make



