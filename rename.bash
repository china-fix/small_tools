#!/bin/bash

rename 's/_1\.fq\.gz/_R1.fq.gz/' *_1.fq.gz
rename 's/_2\.fq\.gz/_R2.fq.gz/' *_2.fq.gz

rename 's/_R1_001\.fastq\.gz/_R1.fq.gz/' *_R1_001.fastq.gz
rename 's/_R2_001\.fastq\.gz/_R2.fq.gz/' *_R2_001.fastq.gz


rename 's/_R1\.fastq\.gz/_R1.fq.gz/' *_R1.fastq.gz
rename 's/_R2\.fastq\.gz/_R2.fq.gz/' *_R2.fastq.gz


rename 's/_1\.fastq\.gz/_R1.fq.gz/' *_1.fastq.gz
rename 's/_2\.fastq\.gz/_R2.fq.gz/' *_2.fastq.gz
