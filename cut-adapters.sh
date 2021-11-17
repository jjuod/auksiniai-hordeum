set -e

# 1. trim adapters

DIR=~/Documents/mieziai/

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o ${DIR}/QCed/ES_11_frw.fastq.gz \
    -p ${DIR}/QCed/ES_11_rev.fastq.gz \
    -m 30 \
    -j 3 \
    ${DIR}/pool_8_11_ES/pool_8_11_ES_UKKD18120008-A12-A44_HF7VKDSXX_L1_1.fq.gz \
    ${DIR}/pool_8_11_ES/pool_8_11_ES_UKKD18120008-A12-A44_HF7VKDSXX_L1_2.fq.gz &> ${DIR}/QCed/ES_11_cutadapt_log.txt
