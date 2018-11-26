#!/usr/bin/env python3


class FileType:
    def __init__(self, name: str, input_variable_name: str, output_variable_name: str):
        self.name: str = name
        self.input_variable_name: str = input_variable_name
        self.output_variable_name: str = output_variable_name


FASTQ = FileType(name="fastq", input_variable_name="InputReads", output_variable_name="OutputReads")
BAM = FileType(name="bam", input_variable_name="InputBams", output_variable_name="OutputBams")
BAI = FileType(name="bai", input_variable_name="InputBais", output_variable_name="OutputBais")
VCF = FileType(name="vcf", input_variable_name="InputVcf", output_variable_name="OutputVcf")
VCFIDX = FileType(name="vcfidx", input_variable_name="InputVcfIdx", output_variable_name="OutputVcfIdx")
RECAL_TABLE = FileType(name="recal_table", input_variable_name="RecalTable", output_variable_name="RecalTable")

