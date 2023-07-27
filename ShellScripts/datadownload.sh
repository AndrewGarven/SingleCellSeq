#!/usr/bin/env bash

#!/bin/bash
#SBATCH --account=def-dmberman
#SBATCH --job-name=sc-datadownload
#SBATCH --qos=privileged
#SBATCH --nodes=2                # number of Nodes
#SBATCH --tasks-per-node=4        # number of MPI processes per node
#SBATCH --mem 8g
#SBATCH --time 00:10:00
#SBATCH --output=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/datadownload.%J.out
#SBATCH --error=/home/garvena/projects/def-dmberman/garvena/singlecellseq/SlurmOut/datadownload.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#

#!/bin/bash

test=(ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_1.fastq.gz
        ftp://ftp.sra.ebi.asqc.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_2.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_2.fastq.gz)

data=(
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/064/SRR17259464/SRR17259464_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/064/SRR17259464/SRR17259464_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/063/SRR17259463/SRR17259463_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/063/SRR17259463/SRR17259463_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/062/SRR17259462/SRR17259462_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/062/SRR17259462/SRR17259462_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/001/SRR9897621/SRR9897621_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/001/SRR9897621/SRR9897621_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/002/SRR9897622/SRR9897622_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/002/SRR9897622/SRR9897622_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/062/SRR12539462/SRR12539462_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/062/SRR12539462/SRR12539462_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/058/SRR14615558/SRR14615558_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/058/SRR14615558/SRR14615558_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/004/SRR9897624/SRR9897624_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/004/SRR9897624/SRR9897624_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/063/SRR12539463/SRR12539463_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/063/SRR12539463/SRR12539463_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/003/SRR9897623/SRR9897623_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/003/SRR9897623/SRR9897623_2.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/005/SRR9897625/SRR9897625_1.fastq.gz"
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/005/SRR9897625/SRR9897625_2.fastq.gz"
)

# Function to download a single URL and save the file with the corresponding file name
parallel --jobs 2 curl -O {} ::: "${test[@]}"

#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_1.fastq.gz -o SRR12603789_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/089/SRR12603789/SRR12603789_2.fastq.gz -o SRR12603789_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_1.fastq.gz -o SRR12603787_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/087/SRR12603787/SRR12603787_2.fastq.gz -o SRR12603787_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_1.fastq.gz -o SRR12603786_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/086/SRR12603786/SRR12603786_2.fastq.gz -o SRR12603786_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_1.fastq.gz -o SRR12603790_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.asqc.uk/vol1/fastq/SRR126/090/SRR12603790/SRR12603790_2.fastq.gz -o SRR12603790_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_1.fastq.gz -o SRR12603788_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/088/SRR12603788/SRR12603788_2.fastq.gz -o SRR12603788_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_1.fastq.gz -o SRR12603784_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/084/SRR12603784/SRR12603784_2.fastq.gz -o SRR12603784_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_1.fastq.gz -o SRR12603785_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/085/SRR12603785/SRR12603785_2.fastq.gz -o SRR12603785_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_1.fastq.gz -o SRR12603783_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/083/SRR12603783/SRR12603783_2.fastq.gz -o SRR12603783_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_1.fastq.gz -o SRR12603782_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/082/SRR12603782/SRR12603782_2.fastq.gz -o SRR12603782_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_1.fastq.gz -o SRR12603781_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/081/SRR12603781/SRR12603781_2.fastq.gz -o SRR12603781_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_1.fastq.gz -o SRR12603780_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/080/SRR12603780/SRR12603780_2.fastq.gz -o SRR12603780_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/064/SRR17259464/SRR17259464_1.fastq.gz -o SRR17259464_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/064/SRR17259464/SRR17259464_2.fastq.gz -o SRR17259464_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/063/SRR17259463/SRR17259463_1.fastq.gz -o SRR17259463_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/063/SRR17259463/SRR17259463_2.fastq.gz -o SRR17259463_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/062/SRR17259462/SRR17259462_1.fastq.gz -o SRR17259462_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/062/SRR17259462/SRR17259462_2.fastq.gz -o SRR17259462_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/001/SRR9897621/SRR9897621_1.fastq.gz -o SRR9897621_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/001/SRR9897621/SRR9897621_2.fastq.gz -o SRR9897621_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/002/SRR9897622/SRR9897622_1.fastq.gz -o SRR9897622_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/002/SRR9897622/SRR9897622_2.fastq.gz -o SRR9897622_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/062/SRR12539462/SRR12539462_1.fastq.gz -o SRR12539462_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/062/SRR12539462/SRR12539462_2.fastq.gz -o SRR12539462_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/058/SRR14615558/SRR14615558_1.fastq.gz -o SRR14615558_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR146/058/SRR14615558/SRR14615558_2.fastq.gz -o SRR14615558_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/004/SRR9897624/SRR9897624_1.fastq.gz -o SRR9897624_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/004/SRR9897624/SRR9897624_2.fastq.gz -o SRR9897624_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/063/SRR12539463/SRR12539463_1.fastq.gz -o SRR12539463_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/063/SRR12539463/SRR12539463_2.fastq.gz -o SRR12539463_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/003/SRR9897623/SRR9897623_1.fastq.gz -o SRR9897623_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/003/SRR9897623/SRR9897623_2.fastq.gz -o SRR9897623_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/005/SRR9897625/SRR9897625_1.fastq.gz -o SRR9897625_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR989/005/SRR9897625/SRR9897625_2.fastq.gz -o SRR9897625_2.fastq.gz
