#!/bin/sh

HOME_BIN=~/bin
# Create dir if $HOME_BIN does not exist.
mkdir -p $HOME_BIN
echo "export PATH=$HOME_BIN:$PATH" >> ~/.bashrc
source ~/.bashrc


echo "This will modify your system installing packages if they are not installed"
read -r -p "Are you sure to install python-biopython cpanminus ncoils? [y/N] " response
case $response in
    [yY][eE][sS]|[yY])
        # install necessary *.deb packages
        sudo apt-get install python-biopython cpanminus ncoils clustalo
        ;;
    *)
        echo "deb packages will not be installed"
        ;;
esac

read -r -p "Are you sure to install Moose Perl module? [y/N] " response
case $response in
    [yY][eE][sS]|[yY])
        # install necessary Perl modules
        cpanm Moose
        ;;
    *)
        echo "Perl module will not be installed"
        ;;
esac

read -r -p "Are you sure to install hmmer with easel miniapps? [y/N] " response
case $response in
    [yY][eE][sS]|[yY])
        # download and compile hmmer with easel miniapps
        cd /tmp
        wget http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
        tar xvf hmmer-3.1b2.tar.gz
        cd hmmer-3.1b2
        ./configure
        make && sudo make install
        cd easel/miniapps
        sudo make install
        ;;
    *)
        echo "hmmer with easel miniapps will not be installed"
        ;;
esac

read -r -p "Are you sure to install Prosite data and executables? [y/N] " response
case $response in
    [yY][eE][sS]|[yY])
		# install Uniprot data and executables
		wget ftp://ftp.expasy.org/databases/prosite/prosite.dat
		mv prosite.dat $HOME_BIN/data
		wget ftp://ftp.expasy.org/databases/prosite/ps_scan/ps_scan_linux_x86_elf.tar.gz
		tar xvf ps_scan_linux_x86_elf.tar.gz
		cd ps_scan
		mv ps_scan.pl pfscan pfsearch psa2msa $HOME_BIN
        ;;
    *)
        echo "Prosite data and executables will not be installed"
        ;;
esac

echo "Assuming you use Bash, this will update your .bashrc"
read -r -p "Are you sure to install PfamScan Perl module? [y/N] " response
case $response in
    [yY][eE][sS]|[yY])
		# install PfamScan Perl module
		cd /tmp
		wget http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz
		tar xvf PfamScan.tar.gz
		mkdir -p ~/perl5/lib/perl5; mv PfamScan $_
		echo "export PERL5LIB=~/perl5/lib/perl5:$PERL5LIB" >> ~/.bashrc
		source .bashrc
        ;;
    *)
        echo "PfamScan will not be installed"
        ;;
esac

echo "Pressing sample Pfam-A database"
hmmpress data/Pfam-A.hmm

echo "You have to download full Pfam-A database and unpack it to "$HOME_BIN"/data"