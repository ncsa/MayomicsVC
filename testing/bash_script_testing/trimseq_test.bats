#!/usr/bin/inputs/env bats
load testing_call

# function teardown() {
#   if [ "$(ls -A outputs)" ]
#   then
#     rm outputs/*
#   fi
#   if [ "$(ls -A test_files)" ]
#   then
#     rm test_files/*
#   fi
#
#   if [ "$(ls -A *.fastq.gz)" ]
#   then
#     rm *.fastq.gz
#   fi
#
#   if [ "$(ls -A *.bam*)" ]
#   then
#     rm *.bam*
#   fi
# }

@test "make sure cutadapt is installed" {
  skip "already tested"
  run apt list cutadapt
  [[ "$output" =~ "[installed]" ]]
}

@test "Test run with no arugment" {
  # skip "already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh
  [ "$status" -eq 1 ]
  [[ "$output" =~ "command line input: `\n`" ]]
  [[ "$output" =~ "No arguments passed." ]]
}

@test "Test help tag" {
  # skip "already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh -h
  echo "$output" > test_files/help_output.txt
  sed -i -e 1,4d test_files/help_output.txt
  test=`diff Verification_files/desired_help_output.txt test_files/help_output.txt`
  [ -z "$test" ]
}

@test "Test with nonexistent option" {
  # skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh -Q garbage
  [[ "$output" =~ "command line input: -Q garbage" ]]
  [[ "$output" =~ "Invalid option: -Q" ]]
}

@test "Test successful run paired end and single end" {
  # skip "Already tested"
  if [[ ! -z `ls outputs` ]]; then
    rm outputs/*
  fi

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  sed -i -e 1,5d outputs/paired_ended_output.cutadapt.log
  test1=`diff Verification_files/desired_paired_ended_output.cutadapt.log outputs/paired_ended_output.cutadapt.log`
  [ -z "$test1" ]
  [[ "$output" =~ "Finished trimming adapter sequences." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/single_end_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r null -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P false -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  sed -i -e 1,5d outputs/single_end_output.cutadapt.log
  test1=`diff Verification_files/desired_single_end_output.cutadapt.log outputs/single_end_output.cutadapt.log`
  [ -z "$test1" ]
  [[ "$output" =~ "Finished trimming adapter sequences." ]]

  [[ `ls outputs` =~ "paired_ended_output.cutadapt.log" ]]
  [[ `ls outputs` =~ "paired_ended_output.trimming.TBD.log" ]]
  [[ `ls outputs` =~ "single_end_output.cutadapt.log" ]]
  [[ `ls outputs` =~ "single_end_output.trimming.TBD.log" ]]
}

@test "Test faulty read 1 options" {
  skip "Already tested"
  tests=(dummy_test_blank.fq dummy_test_text.fq dummy_test_text_with_at.fq)
  output_test_text=("REASON=Input read 1 file garbage_test_files/dummy_test_blank.fq is empty or does not exist." "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'" "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'")
  i=0
  for case in ${tests[@]}
  do
    run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
    -s outputs/test_output -l garbage_test_files/${case} -r \
    ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
    -A ${Input_File_Path}/TruSeqAdaptors.fasta -C \
    /usr/bin/ -t 0 -P true -e inputs/env.txt \
    -F ${Mayomics_path}/src/shell/shared_functions.sh
    if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
    then
        [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
    else
        echo "$output"
        [[ "$output" =~ "${output_test_text[$i]}" ]]
    fi
    let "i=i+1"
  done
}

@test "Test faulty read 2 options" {
  skip "already tested"
  tests=(dummy_test_blank.fq dummy_test_text.fq dummy_test_text_with_at.fq)
  output_test_text=("REASON=Input read 2 file garbage_test_files/dummy_test_blank.fq is empty or does not exist." "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'" "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'")
  i=0
  for case in ${tests[@]}
  do
    run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
    -s outputs/test_output \
    -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
    -r garbage_test_files/${case}  \
    -A ${Input_File_Path}/TruSeqAdaptors.fasta -C \
    /usr/bin/ -t 0 -P true -e inputs/env.txt \
    -F ${Mayomics_path}/src/shell/shared_functions.sh
    if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
    then
        [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
    else
        echo "$output"
        [[ "$output" =~ "${output_test_text[$i]}" ]]
    fi
    let "i=i+1"
  done
}

@test "Test faulty adapter options" {
  skip "Already tested"
  tests=(dummy_test_blank.fq dummy_test_text.fq dummy_test_text_with_gt.fq)
  output_test_text=("REASON=Adapters fasta file garbage_test_files/dummy_test_blank.fq is empty or does not exist." " At line 1: Expected '>' at beginning of FASTA record, but got 'Lorem ipsum dolor sit amet, consectetur adipiscing elit." "is not a valid IUPAC code. Use only characters XACGTURYSWKMBDHVN.")
  i=0
  for case in ${tests[@]}
  do
    run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
    -s outputs/test_output \
    -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
    -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
    -A  garbage_test_files/${case} -C \
    /usr/bin/ -t 0 -P true -e inputs/env.txt \
    -F ${Mayomics_path}/src/shell/shared_functions.sh
    if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
    then
        [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
    else
        [[ "$output" =~ "${output_test_text[$i]}" ]]
    fi
    let "i=i+1"
  done
}

@test "Test bad cutadapt path options" {
  skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/win/ -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Cutadapt directory /usr/win/ is not a directory or does not exist." ]]
}

@test "Test thread value options" {
  skip "Already tested"
  values=(0 32 321 3299 12322)
  for value in ${values[@]}
  do
    run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
    -s outputs/test_output \
    -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
    -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
    -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
    /usr/bin/ -t ${value} -P true -e inputs/env.txt \
    -F ${Mayomics_path}/src/shell/shared_functions.sh
    if [ $value -lt 321 ]
    then
      [[ "$output" =~ "Finished trimming adapter sequences." ]]
    else
      [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
    fi
  done
}

@test "Test true/false values" {
  skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 12 -P True -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 12 -P T -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r null \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 12 -P False -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r null \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 12 -P f -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]
}

@test "Incorrect read options" {
  skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -A ${Input_File_Path}/TruSeqAdaptors.fasta \
  -C /usr/bin/ -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Missing read 2 option: -r. If running a single-end job, set -r null in command." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
   -A ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 0 -P false -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=User specified Single End option, but did not set read 2 option -r to null." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r null \
  -A ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Input read 2 file null is empty or does not exist." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/test_output \
  -l null -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -A ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  /usr/bin/ -t 0 -P false -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Input read 1 file null is empty or does not exist." ]]

}

@test "Test missing option values" {
  skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -s in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -l in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -r in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -A in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -C in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -t in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -P in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Error with option -e in command. Option passed incorrectly or without argument." ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F
  [[ "$output" =~ "Option -F requires an argument." ]]
}

@test "Test path name errors" {
  skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq. \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Input read 1 file" ]]
  [[ "$output" =~ "WGS_chr1_5X_E0.005_L1_read1.fastq. is empty or does not exist" ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq. \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Input read 2 file" ]]
  [[ "$output" =~ "WGS_chr1_5X_E0.005_L1_read2.fastq. is empty or does not exist" ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fsta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Adapters fasta file" ]]
  [[ "$output" =~ "TruSeqAdaptors.fsta is empty or does not exist" ]]

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inpats/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "inpats/env.txt: No such file or directory" ]]
}

@test "Checks file permissions" {
  skip "Already tested"
  run chmod 000 ${Input_File_Path}/inputs
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "REASON=Input read 1 file" ]]
  [[ "$output" =~ "is empty or does not exist." ]]
  run chmod 755 ${Input_File_Path}/inputs

  run chmod 000 outputs
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Permission denied" ]]
  run chmod 755 outputs

  run chmod 000 ${Mayomics_path}/src
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ "$output" =~ "Permission denied" ]]
  run chmod 755 ${Mayomics_path}/src

}

@test "Checks that output file/folder permissions are set" {
  skip "Already tested"
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ `ls -l outputs/paired_ended_output.cutadapt.log` =~ "-rw-r--r--" ]]
# Original
#  [[ `ls -l outputs/paired_ended_output.read1.trimmed.fq.gz` =~ "-rw-r--r--" ]]
  [[ `ls -l WGS_chr1_5X_E0.005_L1_read1.fastq.gz` =~ "-rw-r--r--" ]]
# Original
#  [[ `ls -l outputs/paired_ended_output.read2.trimmed.fq.gz` =~ "-rw-r--r--" ]]
  [[ `ls -l WGS_chr1_5X_E0.005_L1_read2.fastq.gz` =~ "-rw-r--r--" ]]

  run chmod 200 WGS_chr1_5X_E0.005_L1_read1.fastq.gz
  echo `ls -l outputs` >> output.txt
  run //bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
# Original
  # [[ `ls -l outputs/paired_ended_output.read1.trimmed.fq.gz` =~ "--w-r-----" ]]
  [[ `ls -l WGS_chr1_5X_E0.005_L1_read1.fastq.gz` =~ "--w-r-----" ]]
  run chmod 644 WGS_chr1_5X_E0.005_L1_read1.fastq.gz

  run chmod 200 WGS_chr1_5X_E0.005_L1_read1.fastq.gz
  echo `ls -l outputs` >> output.txt
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ `ls -l WGS_chr1_5X_E0.005_L1_read1.fastq.gz` =~ "--w-r-----" ]]
  run chmod 644 WGS_chr1_5X_E0.005_L1_read1.fastq.gz

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  run chmod 200 WGS_chr1_5X_E0.005_L1_read1.fastq.gz

  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/paired_ended_output \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A  ${Input_File_Path}/TruSeqAdaptors.fasta -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh
  [[ `ls -l WGS_chr1_5X_E0.005_L1_read1.fastq.gz` =~ "--w-r-----" ]]
  run chmod 644 WGS_chr1_5X_E0.005_L1_read1.fastq.gz
}

@test "Logs are truncated at the beginning of the run" {
  skip "Already tested"
  # Create a log with a failed run
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/first_run \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A garbage_test_files/dummy_test_text.fq -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh

  IFS=$'\r\n' GLOBIGNORE='*' command eval  'RUN1_CUTLOG=($(cat outputs/first_run.cutadapt.log))'
  sed -i -e 1,1d outputs/first_run.trimming.TBD.log
  sed -i -e 4,4d outputs/first_run.trimming.TBD.log
  IFS=$'\r\n' GLOBIGNORE='*' command eval  'RUN1_STDOUT=($(cat outputs/first_run.trimming.TBD.log))'

  # Second run has a slightly different problem
  run /bin/bash ${Mayomics_path}/src/shell/trim_sequences.sh \
  -s outputs/second_run \
  -l ${Input_File_Path}/inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz \
  -r ${Input_File_Path}/WGS_chr1_5X_E0.005_L1_read2.fastq.gz \
  -A garbage_test_files/dummy_test_text_with_gt.fq -C /usr/bin/ \
  -t 0 -P true -e inputs/env.txt \
  -F ${Mayomics_path}/src/shell/shared_functions.sh


  # Check to make sure that the logs are overwritten and don't contain old information
  sed -i -e 1,1d outputs/second_run.trimming.TBD.log
  sed -i -e 4,4d outputs/second_run.trimming.TBD.log
  [[ `cat outputs/second_run.cutadapt.log` =~ "cutadapt: error: Character 'I' in adapter sequence 'VIVAMTS TLTRICES FELIS AT METTS MALESTADA IMPERDIET.STSPENDISSE A LEO BLANDIT, CONSEQTAT EX ET, FATCIBTS TTRPIS.ETIAM MOLLIS RISTS QTIS ERAT ELEIFEND BIBENDTM.STSPENDISSE SODALES MAGNA ID LIGTLA SAGITTIS, SIT AMET LOBORTIS NIBH POSTERE.PHASELLTS VENENATIS NEQTE AT NEQTE PELLENTESQTE, ET ELEMENTTM MI TLTRICIES.' is not a valid IUPAC code. Use only characters XACGTURYSWKMBDHVN." ]]
  # echo "${RUN1_CUTLOG[13]}" >&3
  [[ ! `cat outputs/second_run.cutadapt.log` =~ "${RUN1_CUTLOG[13]}" ]]
  [[ `cat outputs/second_run.trimming.TBD.log` =~ "trim_sequences.sh stopped at line 260. Cutadapt Read 1 and 2 failure." ]]
  [[ ! `cat outputs/second_run.trimming.TBD.log` =~ "${RUN1_STDOUT[3]}" ]]

}

# @test ""
