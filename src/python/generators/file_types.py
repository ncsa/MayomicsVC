#!/usr/bin/env python3

from typing import List


class FileType:
    def __init__(self, name: str, extensions: List[str]):
        self.name: str = name
        self.extensions: List[str] = extensions
        self.inUse: bool = False

    def markAsUsed(self):
        self.inUse = True


FASTQ = FileType(name="fastq", extensions=[".fq", ".fastq"])
BAM = FileType(name="bam", extensions=[".bam"])
BAI = FileType(name="bai", extensions=[".bai"])
VCF = FileType(name="vcf", extensions=[".vcf"])
RECAL_TABLE = FileType(name="recal_table", extensions=[".table"])
