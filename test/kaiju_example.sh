#!/bin/sh

#### STEP 0: Setup paths
#
KAIJU_HOME=$HOME/proj/SDK/sdk_modules/kb_kaiju/bin/kaiju
KAIJU_BINDIR=$KAIJU_HOME/bin
KRONA_HOME=$HOME/proj/SDK/sdk_modules/kb_kaiju/bin/Krona
KRONA_BINDIR=$KRONA_HOME/bin/bin

KAIJU_BIN=$KAIJU_BINDIR/kaiju
KAIJU_REPORT_BIN=$KAIJU_BINDIR/kaiju2table
KAIJU2KRONA_BIN=$KAIJU_BINDIR/kaiju2krona
KRONAIMPORT_BIN=$KRONA_BINDIR/ktImportText

KAIJU_DBDIR=$HOME/proj/SDK/sdk_modules/kb_kaiju/data/kaijudb
KAIJU_DBTYPE=refseq
#KAIJU_DBTYPE=progenomes
#KAIJU_DBTYPE=nr
#KAIJU_DBTYPE=nr_euk
KAIJU_NODES=$KAIJU_DBDIR/$KAIJU_DBTYPE/nodes.dmp
KAIJU_NAMES=$KAIJU_DBDIR/$KAIJU_DBTYPE/names.dmp

if [ $KAIJU_DBTYPE = "refseq" ] ; then
    KAIJU_DBPATH=$KAIJU_DBDIR/$KAIJU_DBTYPE/kaiju_db_refseq.fmi
elif [ $KAIJU_DBTYPE = "progenomes" ] ; then
    KAIJU_DBPATH=$KAIJU_DBDIR/$KAIJU_DBTYPE/kaiju_db_progenomes.fmi
elif [ $KAIJU_DBTYPE = "nr" ] ; then
    KAIJU_DBPATH=$KAIJU_DBDIR/$KAIJU_DBTYPE/kaiju_db_nr.fmi
elif [ $KAIJU_DBTYPE = "nr_euk" ] ; then
    KAIJU_DBPATH=$KAIJU_DBDIR/$KAIJU_DBTYPE/kaiju_db_nr_euk.fmi
fi

RUN_DATA_DIR=$HOME/proj/SDK/sdk_modules/kb_kaiju/test/data
fwd_reads=$RUN_DATA_DIR/seven_species_nonuniform_10K.PE_reads_fwd-0.FASTQ.gz
#fwd_reads=$RUN_DATA_DIR/seven_species_nonuniform_100K.PE_reads_fwd-0.FASTQ.gz
rev_reads=$RUN_DATA_DIR/seven_species_nonuniform_10K.PE_reads_rev-0.FASTQ.gz
#rev_reads=$RUN_DATA_DIR/seven_species_nonuniform_100K.PE_reads_rev-0.FASTQ.gz
out_basename=test2_multi_tax
kaiju_out_file=$RUN_DATA_DIR/$out_basename.kaiju.out
kaiju_summary_out_file=$kaiju_out_file.summary
krona_in_file=$kaiju_out_file.krona
krona_out_html_file=$kaiju_out_file.krona.html


#### STEP 1: kaiju run and summary
#
threads="-z 4"
minlength="-m 11"
minscore="-s 65"
mismatches="-e 5"
e_value="-E 0.05"
verbose="-v"
greedy="-a greedy"
SEG_filter="-x"
fwd_reads_arg="-i $fwd_reads"
if [ -s rev_reads ] ; then
    rev_reads_arg="-j $rev_reads"
else
    rev_reads_arg=""
fi
cmd="$KAIJU_BIN -t $KAIJU_NODES -f $KAIJU_DBPATH -i $fwd_reads_arg -j $rev_reads_arg -o $kaiju_out_file $SEG_filter $minlength $minscore $greedy $mismatches $e_value $threads $verbose"
if [ ! -s $kaiju_out_file ] ; then
    echo $cmd
    exec $cmd
fi

# kaiju summary run
#tax_levels="phylum,class,order,family,genus,species"
tax_level="family"
#taxon_fullpath_arg="-p"
taxon_fullpath_arg=""
filter_perc=1
filter_unclassified="-u"
#filter_unclassified=""
if [ $filter_perc -gt 0 ] ; then
    filter_arg="-m $filter_perc $filter_unclassified"
else
    filter_arg=$filter_unclassified
fi
cmd="$KAIJU_REPORT_BIN -t $KAIJU_NODES -n $KAIJU_NAMES -r $tax_level $filter_arg $taxon_fullpath_arg -o $kaiju_summary_out_file $kaiju_out_file"
if [ ! -s $kaiju_summary_out_file ] ; then
    echo
    echo $cmd
    exec $cmd
fi


#### STEP 2: krona prep and html gen run
#
cmd="$KAIJU2KRONA_BIN -t $KAIJU_NODES -n $KAIJU_NAMES -i $kaiju_out_file -o $krona_in_file"
if [ ! -s $krona_in_file ] ; then
    echo
    echo $cmd
    exec $cmd
fi

cmd="$KRONAIMPORT_BIN -o $krona_out_html_file $krona_in_file"
if [ ! -s $krona_out_html_file ] ; then
    echo
    echo $cmd
    exec $cmd
fi


#### STEP 3: Finalize
#
echo "DONE"
exit 0
