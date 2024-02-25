#!/bin/sh

#set -x

module load jdk/11.0.2
#module load nextflow_latest/20.10.0
NEXTFLOW_VER=$(realpath ../../nextflow_ver/20.10.0)
NEXTFLOW_DIR=$(realpath ../../nextflow)
PARAMS_FILE= 
CONFIGURATION_FILE=
OUTPUT_DIR=
EXECUTOR= 
PARAMS_BASENAME= 
DEMUX_DIR=
CONFIG_DIR=
SEGMENT="12"
CURRENT_DIR=$PWD

usage() {
        echo "
Program: Wrapper(CLI) for the bfx-ngs pipeline
Usage:
        bash run_nextflow_demux_single-sample.sh [options]
Options: -p|--params         parameter file
         -o|--output-dir     output directory
         -e|--executor       executor (local|cluster|aws)
         -seg|--segment      semgment (demux: demux only | single-sample: sample process only). (Default=12)
         -d|--demux-dir      demultiplexed directory only if -seg=single-sample. (Default=$OUTPUT_DIR/analysis/demux)
         -h|--help
"
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -p|--params)
            PARAMS_FILE="$(realpath $2)"
            PARAMS_BASENAME=`basename $PARAMS_FILE`
            shift
        ;;
        -o|--output-dir)
            OUTPUT_DIR="$(realpath $2)"
            DEMUX_DIR=$OUTPUT_DIR/analysis/demux
            CONFIG_DIR=$OUTPUT_DIR/config
            shift
        ;;
        -e|--executor)
            EXECUTOR="$2"
            shift
        ;;
        -seg|--segment)
            SEGMENT="$2"
            shift
        ;;
        -d|--demux-dir)
            DEMUX_DIR="$(realpath $2)"
            shift
        ;;
        -h|--help)
            usage
            exit 1
        ;;
     esac
     shift
done

echo "PARAMS_FILE : $PARAMS_FILE"
echo "CONFIGURATION_FILE : $CONFIGURATION_FILE"
echo "OUTPUT_DIR : $OUTPUT_DIR"
echo "DEMUX_DIR : $DEMUX_DIR"

if [[ $SEGMENT == 'single-sample' ]] && [[ ! -d $DEMUX_DIR ]]; then
    (>&2 echo "demuxed directory does not exists or invalid path.")
    usage
    exit 1
fi

demux() {
    echo "running demultiplexing ..."
    cd $CONFIG_DIR/demux

    cmd="$NEXTFLOW_VER/nextflow run $NEXTFLOW_DIR/workflows/demultiplex.nf -ansi-log false -profile $EXECUTOR -c $PARAMS_BASENAME | tee -a log.txt"
    echo "$cmd" > run_nextflow.sh && eval "$cmd"
    cd $CURRENT_DIR
}

setup_demux(){
    mkdir -p $CONFIG_DIR/demux
    cp $PARAMS_FILE $CONFIG_DIR/demux/
    echo "/// params for demux" >> $CONFIG_DIR/demux/$PARAMS_BASENAME
    echo "params.finalDir = \""$OUTPUT_DIR"\"" >> $CONFIG_DIR/demux/$PARAMS_BASENAME
}

single_sample() {
    local sample_name=$1
    echo "running $sample_name ..."
    cd $CONFIG_DIR/$sample_name

    cmd="$NEXTFLOW_VER/nextflow run $NEXTFLOW_DIR/workflows/single-sample.nf -ansi-log false -profile $EXECUTOR -c $PARAMS_BASENAME -w $CONFIG_DIR/$sample_name/work | tee -a log.txt"
    echo -e "export NXF_ANSI_LOG=false\n$cmd" >> run_nextflow.sh
    eval "$cmd"
    cd $CURRENT_DIR
}

setup_single_sample() {
    local sample_name=$1

    mkdir -p $CONFIG_DIR/$sample_name
    cp $PARAMS_FILE $CONFIG_DIR/$sample_name/
    echo "/// params for single_sample" >> $CONFIG_DIR/$sample_name/$PARAMS_BASENAME
    echo "params.sampleName = \""$sample_name"\"" >> $CONFIG_DIR/$sample_name/$PARAMS_BASENAME
    echo "params.sampleFastq = \""$DEMUX_DIR/$sample_name/$(ls $DEMUX_DIR/$sample_name)"\"" >> $CONFIG_DIR/$sample_name/$PARAMS_BASENAME 
    echo "params.finalDir = \""$OUTPUT_DIR"\"" >> $CONFIG_DIR/$sample_name/$PARAMS_BASENAME
}

##
##- the main routine
##

if [[ ! $SEGMENT == 'single-sample' ]]; then
    setup_demux && demux || exit $?
fi

if [[ ! $SEGMENT == 'demux' ]]; then
    for f in `find $DEMUX_DIR -mindepth 1 -type d -exec basename {} \;`; 
    do 
        setup_single_sample $f || exit $?
        single_sample $f &
        pids[${f}]=$!
        sleep 1
    done

    # wait for all pids
    for pid in ${pids[*]}; do
        wait $pid
    done
fi

