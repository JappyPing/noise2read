#!/bin/bash
# @Author: Pengyao Ping
# @Date:   2022-12-19 18:51:21
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-27 23:24:17
#!/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
core_count=$(($( grep -c ^processor /proc/cpuinfo)-2))
echo "using $core_count threads"


#set -e
#set -o errexit

#---extract commandline arguments---
function  help {
    echo "Usage:"
    echo "Typical use for paired-end reads, align both read files (-1 read1 -2 read2) together (-m tog) with \$aligner (-a \$aligner):"
    echo "run.sh -r reference.fasta -1 first_read_file.fastq -2 second_read_file.fastq"
    echo "-o output_directory"
    exit 2
}

out_dir="none"

while getopts ":h:r:1:2:o:" option; do
    ((optnum++))
   case ${option} in
      h) # display Help

         help
         exit;;
     \?)
         echo "Invalid option, use run.sh -h for usage"
         exit;;
     r)
         ref=${OPTARG};;
     1)
         read1=${OPTARG};;
     2)
         read2=${OPTARG};;
    o)
        out_dir=${OPTARG};;

   esac
done


if [ $OPTIND -eq 1 ]; then
    help
fi
shift $((OPTIND -1))

#--- build index----
function build_index {
    if ! test -f $ref;then
        echo "reference file $ref not existed"
        exit
    fi
    echo "building index for $ref"

    if [ -d ${DIR}/bowtie2_index ];then
        rm ${DIR}/bowtie2_index/*
    else
        mkdir ${DIR}/bowtie2_index
    fi
    ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2-build --large-index $ref ${DIR}/bowtie2_index/reference_index

}

#---alignment---

function  alignment {

    if [ -e "$read1" ] && [ -e "$read2" ];then
        echo "bowtie2 alignment"
        ${DIR}/bowtie2-2.4.4-linux-x86_64/bowtie2 -p ${core_count} -x ${DIR}/bowtie2_index/reference_index -1 "$read1" -2 "$read2" --local > gene.sam

    else
        echo "Error getting read files, use -0 for single read file, use -1 read1.fq -2 read2.fq for paired-end reads"
        exit
    fi
}

function split_sam {
    gene_size=$[$(stat --printf="%s" $1)]
    if [ "$gene_size" -eq 0 ]; then
        echo "error empty alignment result, exiting"
        exit
    fi

    #--- sample about 600 reads to estimate the lines of gene.sam useful for very large gene.sam ---
    sample_size=300
    head -$((sample_size+2)) $1|tail -$sample_size > temp_size.sam && tail -$sample_size $1 >> temp_size.sam
    avg=$[(($(ls -l temp_size.sam | awk '{print $5}')/(sample_size*2)))]
    estimate=$[(($(ls -l $1 | awk '{print $5}')/avg))]
    divide=$((estimate/core_count))

    echo $divide $estimate $avg
    if [ -z $divide ] || [ -z $estimate ] || [ -z $avg ]; then
        echo "error values, d $divide e $estimate a $avg c $core_count"
        exit
    fi
    time split -l $divide $1 split_gene

    rm temp_size.sam

}

function extract_sam {

    touch extract.sam
    cp /dev/null extract.sam
    cp /dev/null full_extract.sam
    count=0
    arr=()

    for f in split_gene*
    do
        touch extract_$count.sam
        cp /dev/null extract_$count.sam
        cp /dev/null full_extract_$count.sam
        echo "split gene.sam into files $count"
        awk 'BEGIN{tmp="'"*"'"}{if(NF>10 && $6!=tmp) print $0}' $f  > full_extract_$count.sam  &
        arr+=($!)
        count=$((count+1))
    done
    for ((x=0;x<=$count;x++))
    do
        wait ${arr[$x]}
    done

    arr=()

    x=0
    for full1 in full_extract_*
    do
        cat $full1 >> full_extract.sam
        awk '{ print $1" "$2 " "$4 " "$10" "$6}' $full1  > extract_$x.sam  &
        arr+=($!)
        x=$((x+1))
    done

    for ((x=0;x<=$count;x++))
    do
        echo "wait, $(wc -l "full_extract_$x.sam") lines in full_extract_$x.sam "
        wait ${arr[$x]}
    done
    echo "extract_ ready"
    for f1 in extract_*
    do
        echo "$[$(grep -c '^' $f1 )] lines from $f1 into extract.sam"
        cat $f1 >>extract.sam
    done

    rm split_gene*
    rm extract_*
    rm full_extract_*
    if [   $[$(wc -l extract.sam|cut -d' ' -f1)] -eq 0 ];then
        echo "no matched reads"
        exit
    fi
    cat full_extract.sam | sort -nk4 > s_full_extract.sam
    cp s_full_extract.sam full_extract.sam
    cat extract.sam | sort -nk3 > s_extract.sam
    cat s_extract.sam > extract.sam
    #remove firstline of bowtie2 output that causes error
    sed -i 1d extract.sam
    echo  $[$(wc -l extract.sam|cut -d' ' -f1)] " total reads"
    printf "\n\n produce extract \n rl \"${rL}\" \n"
    echo " extract.sam  with" $[$(wc -l extract.sam|cut -d' ' -f1)] "reads"
}

if  [ "$out_dir" == "none" ]; then
    out_dir="${ref##*/}"_output
fi

echo output dir is  "$out_dir"


fL=$( tr -d '\r' < "$ref" | awk 'BEGIN{RS="\n"}{if(substr($0,1,1) ~ /\>/ ) printf("%s", $NF);next}  '  | wc -m)
if (( $(echo "$fL < $len_limit" | bc -l) )); then
    echo "gene length $fL too short"
    exit
fi

if [ -d "$out_dir" ]; then
    rm -r "$out_dir"
    mkdir "$out_dir"
else
    mkdir "$out_dir"
fi


build_index
alignment

echo "get matched reads"
split_sam gene.sam
extract_sam
rm s_extract.sam gene.sam s_full_extract.sam


python3 "${DIR}"/get_coverage.py "${ref}" extract.sam

mv *.sam ${out_dir}

mv *.txt ${out_dir}
