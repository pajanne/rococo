#!/bin/sh

# --------------------------------------------------------
# Set all the environment variables for tracking scripts
# nds@sanger.ac.uk, ap12@sanger.ac.uk
# --------------------------------------------------------

umask 002

# Where is the tracking code?
export SOFT_VERTRES=/software/vertres
export VERTRES_SVN=/nfs/pathogen/svn/trunk
export VERTRES_BIN=$VERTRES_SVN/bin


# LSF variables
export LSB_DEFAULTPROJECT=team81
export LSF_BINDIR=/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/bin
export LSF_ENVDIR=/etc
export LSF_LIBDIR=/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/lib
export LSF_SERVERDIR=/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/etc
export XLSF_UIDDIR=/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/lib/uid

# mpsa_download
export ORACLE=/software/oracle
export ORACLE_HOME=/software/oracle
export LD_LIBRARY_PATH=/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/lib:/software/badger/lib:/software/oracle/lib

# perl
export PERL_INLINE_DIRECTORY=$SOFT_VERTRES/lib/inline
export PERL5LIB=$SOFT_VERTRES/lib/all:$VERTRES_SVN/modules:$VERTRES_BIN:/nfs/users/nfs_n/nds/nds/:/nfs/users/nfs_n/nds/nds/googledocs/classes
export PERLDOC_PAGER=less

# python
export PYTHONPATH=/nfs/users/nfs_a/ap12/genlibpy:/software/pathogen/external/lib/python/lib/python2.6/site-packages/:/software/pathogen/psu_svn/trunk/genlib/python

# path
export PATH=/software/pathogen/external/lib/python/bin/:/software/pathogen/external/bin/:/software/pathogen/external/applications/EMBOSS/bin/:$ORACLE_HOME/bin:$VERTRES_BIN:/software/solexa/bin:/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/etc:/usr/local/lsf/7.0/linux2.6-glibc2.3-x86_64/bin:/software/bin:/usr/local/bin:/usr/bin:/bin:/nfs/users/nfs_n/nds/nds/googledocs/classes

# tools
export SAMTOOLS=$VERTRES_BIN/samtools-0.1.7
export GATK=$VERTRES_BIN/GenomeAnalysisTK
export GATK_RESOURCES=/lustre/scratch102/g1k/ref/broad_recal_data
export PICARD=$VERTRES_BIN/picard-tools
export DINDEL_SCRIPTS=$VERTRES_BIN/dindel-0.1

# tracking databases
export VRTRACK_HOST=web-mii-shap
export VRTRACK_PORT=3303
export VRTRACK_RO_USER=pathpipe_ro
export VRTRACK_RW_USER=pathpipe_rw
export VRTRACK_PASSWORD=path3476

# inline java, proxy and google document details, currently pointing to java classes in nds's directory (will change later)
export CLASSPATH=/nfs/users/nfs_n/nds/nds/googledocs/lib/gdata-core-1.0.jar:/nfs/users/nfs_n/nds/nds/googledocs/lib/gdata-spreadsheet-3.0.jar:/nfs/users/nfs_n/nds/nds/googledocs/lib/google-collect-1.0-rc1.jar:.:/nfs/users/nfs_n/nds/nds:/nfs/users/nfs_n/nds/nds/googledocs/classes

export HTTP_PROXY="wwwcache.sanger.ac.uk"
export HTTP_PORT="3128"
export HTTPS_PROXY="wwwcache.sanger.ac.uk"
export HTTPS_PORT="3128"

export google_user="pathpipe@gmail.com"
export google_password="pathpipeE204"
export google_doc="path_organism_list"

# data hierarchy structure
export DATA_HIERARCHY="genus:species-subspecies:TRACKING:projectid:sample:technology:library:lane"
