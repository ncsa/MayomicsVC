#!/usr/bin/inputs/env bats
load testing_call
# Instructions: run using full path from within your Mayomics/src directory
# TODO update to make the above true

# @test "make sure cutadapt is installed" {
#     # skip "already tested"
#     run apt list cutadapt
#     [[ "$output" =~ "[installed]" ]]
# }
#
# @test "Test no-arugment, fake argument, and help functions" {
#     # skip "already tested"
#     run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh
#     [ "$status" -eq 1 ]
#     echo "$output" > test_files/noarg_output.txt
#     sed -i -e 1,4d test_files/noarg_output.txt
#     test=`diff Verification_files/desired_noarg_output.txt test_files/noarg_output.txt`
#
#     run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh -h
#     echo "$output" > test_files/help_output.txt
#
#     sed -i -e 1,4d test_files/help_output.txt
#     test=`diff Verification_files/desired_help_output.txt test_files/help_output.txt`
#     [ -z "$test" ]
#
#     run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh -Q garbage
#     echo "$output" > test_files/nonexist_output.txt
#
#     sed -i -e 1,4d test_files/nonexist_output.txt
#     test=`diff Verification_files/desired_nonexist_output.txt test_files/nonexist_output.txt`
#     [ -z "$test" ]
# }
#
# @test "Test faulty read 1 options" {
#     tests=(dummy_test_blank.fq dummy_test_text.fq dummy_test_text_with_at.fq)
#     output_test_text=("REASON=Input read 1 file garbage_test_files/dummy_test_blank.fq is empty or does not exist." "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'" "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'")
#     i=0
#     for case in ${tests[@]}
#     do
#       run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#       -s outputs/test_output -l garbage_test_files/${case} -r \
#       inputs/Chr1_read1.fq -A inputs/toy_adapters.fa -C \
#       /usr/bin/ -t 0 -P true -e inputs/env.txt
#       if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
#       then
#           [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
#       else
#           echo "$output"
#           [[ "$output" =~ "${output_test_text[$i]}" ]]
#       fi
#       let "i=i+1"
#     done
# }
#
# @test "Test faulty read 2 options" {
#     output_test_text=("REASON=Input read 2 file garbage_test_files/dummy_test_blank.fq is empty or does not exist." "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'" "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'")
#     i=0
#     for case in ${tests[@]}
#     do
#       run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#       -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#       garbage_test_files/${case} -A inputs/toy_adapters.fa -C \
#       /usr/bin/ -t 0 -P true -e inputs/env.txt
#       if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
#       then
#           [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
#       else
#           echo "$output"
#           [[ "$output" =~ "${output_test_text[$i]}" ]]
#       fi
#       let "i=i+1"
#     done
# }
#
# @test "Test faulty adapter options" {
#     #case 3: faulty toy adapter
#     output_test_text=("REASON=Adapters fasta file garbage_test_files/dummy_test_blank.fq is empty or does not exist." "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'" "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'")
#     i=0
#     for case in ${tests[@]}
#     do
#       run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#       -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#       inputs/Chr1_read2.fq -A garbage_test_files/${case} -C \
#       /usr/bin/ -t 0 -P true -e inputs/env.txt
#       if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
#       then
#           [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
#       else
#           echo "$output"
#           [[ "$output" =~ "${output_test_text[$i]}" ]]
#       fi
#       let "i=i+1"
#     done
# }
#
# @test "Test bad cutadapt path options" {
#     run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#     -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#     inputs/Chr1_read2.fq -A inputs/toy_adapters.fa -C \
#     /usr/win/ -t 0 -P true -e inputs/env.txt
#     [[ "$output" =~ "REASON=Cutadapt directory /usr/win/ is not a directory or does not exist." ]]
# }
#
# @test "Test thread value options" {
#     values=(0 1 32 321 3299 12322)
#     for value in ${values[@]}
#     do
#       run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#       -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#       inputs/Chr1_read2.fq -A inputs/toy_adapters.fa -C \
#       /usr/bin/ -t $value -P true -e inputs/env.txt
#       if [ $value -lt 321 ]
#       then
#         [[ "$output" =~ "Finished trimming adapter sequences." ]]
#       else
#         [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
#       fi
#     done
# }
#
# @test "Test faulty adapter options" {
#     output_test_text=("REASON=Adapters fasta file garbage_test_files/dummy_test_blank.fq is empty or does not exist." "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'" "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'")
#     i=0
#     for case in ${tests[@]}
#     do
#       run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#       -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#       inputs/Chr1_read2.fq -A garbage_test_files/${case} -C \
#       /usr/bin/ -t 0 -P true -e inputs/env.txt
#       if [[ "$output" =~ "Cutadapt Read 1 and 2 failure." ]]
#       then
#           [[ `grep "${output_test_text[$i]}" outputs/test_output.cutadapt.log` ]]
#       else
#           echo "$output"
#           [[ "$output" =~ "${output_test_text[$i]}" ]]
#       fi
#       let "i=i+1"
#     done
# }
#
# @test "Test true/false values" {
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#   inputs/Chr1_read2.fq -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P True -e inputs/env.txt
#   [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#   inputs/Chr1_read2.fq -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P t -e inputs/env.txt
#   [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#   inputs/Chr1_read2.fq -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P False -e inputs/env.txt
#   [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#   inputs/Chr1_read2.fq -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P f -e inputs/env.txt
#   [[ "$output" =~ "REASON=Incorrect argument for paired-end option -P. Must be set to true or false." ]]
# }
#
# @test "Incorrect right read options" {
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -A inputs/toy_adapters.fa \
#   -C /usr/bin/ -t 0 -P true -e inputs/env.txt
#   [[ "$output" =~ "REASON=Missing read 2 option: -r. If running a single-end job, set -r null in command." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#    -A inputs/toy_adapters.fa \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P false -e inputs/env.txt
#   [[ "$output" =~ "REASON=User specified Single End option, but did not set read 2 option -r to null." ]]
# }
#
# @test "Test missing option values" {
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P true -e inputs/env.txt
#   [[ "$output" =~ "Error with option -s in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l -r inputs/Chr1_read2.fq \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P false -e inputs/env.txt
#   [[ "$output" =~ "Error with option -l in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P false -e inputs/env.txt
#   [[ "$output" =~ "Error with option -r in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#   -A -C \
#   /usr/bin/ -t 0 -P false -e inputs/env.txt
#   [[ "$output" =~ "Error with option -A in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#   -A inputs/toy_adapters.fa -C \
#   -t 0 -P false -e inputs/env.txt
#   [[ "$output" =~ "Error with option -C in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t -P false -e inputs/env.txt
#   [[ "$output" =~ "Error with option -t in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P -e inputs/env.txt
#   [[ "$output" =~ "Error with option -P in command. Option passed incorrectly or without argument." ]]
#
#   run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
#   -s outputs/test_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
#   -A inputs/toy_adapters.fa -C \
#   /usr/bin/ -t 0 -P false -e
#   [[ "$output" =~ "Option -e requires an argument." ]]
#
# }

@test "Test missing option values" {
  run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
  -s outputs/paired_ended_output -l inputs/Chr1_read1.fq -r inputs/Chr1_read2.fq \
  -A inputs/toy_adapters.fa -C /usr/bin/ -t 0 -P true -e inputs/env.txt
  sed -i -e 5,5d outputs/paired_ended_output.cutadapt.log
  sed -i -e 1,4d outputs/paired_ended_output.trimming.TBD.log
  test1=`diff Verification_files/desired_paired_ended_output.cutadapt.log outputs/paired_ended_output.cutadapt.log`
  test2=`diff Verification_files/desired_paired_ended_output.trimming.TBD.log outputs/paired_ended_output.trimming.TBD.log`
  [ -z "$test1" ]
  [ -z "$test2" ]

  # run /bin/bash ${MAYOMICS_PATH}/src/bash_scripts/trim_sequences.sh \
  # -s outputs/single_end_output -l inputs/Chr1_read1.fq -r null \
  # -A inputs/toy_adapters.fa -C /usr/bin/ -t 0 -P false -e inputs/env.txt
  # echo "$output" > test_files/single_end_output.txt
  # sed -i -e 1,4d test_files/single_end_output.txt
  # test1=`diff Verification_files/desired_single_end_output.cutadapt.log test_files/single_end_output.cutadapt.log`
  # test2=`diff Verification_files/desired_single_end_output.trimming.TBD.log test_files/single_end_output.trimming.TBD.log`
  # [ -z "$test1" ]
  # [ -z "$test2" ]
}
