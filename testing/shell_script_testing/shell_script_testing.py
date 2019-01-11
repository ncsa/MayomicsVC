import sys
import os
import subprocess
import unittest
import time
import glob
import copy


class Script:
    """
    Trim_sequences is done, and the other scripts will have similar setups that have the flags as attributes so tests
    can step through the list of attributes.
    """
    __cwd = os.getcwd()

    def __init__(self):
        self.path = self.__cwd
        self.shell_path = '{}/MayomicsVC/src/shell'.format(self.__cwd)
        self.test_path = '{}/MayomicsVC/testing/shell_script_testing'.format(self.__cwd)


class Trimming(Script):
    """
    Constructs the trim_sequences.sh commands for single and paired end reads. Each 'flag' attribute represents a
    particular flag, that way we can step through the flags and perform tests on each.

    TODO: os.path to make tests more generic
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s {}/outputs/output".format(self.path)
        self.flag_A = '-A Inputs/TruSeqAdaptors.fasta'
        self.flag_l = '-l Inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz'
        self.flag_r = '-r Inputs/WGS_chr1_5X_E0.005_L1_read2.fastq.gz'
        self.flag_C = '-C /usr/local/apps/bioapps/python/Python-3.6.1/bin'  # for iforge testing
        # self.flag_C = '-C /usr/bin' # for local testing
        self.flag_t = '-t 8'
        self.flag_P = '-P true'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/trim_sequences.sh'.format(self.shell_path)
        self.type = 'trim_sequences.sh'

    def __str__(self, case: str = 'paired'):
        if case == 'single':
            return "/bin/bash {} {} {} {} -r null {} {} -P false {} {} {}".format(self.name, self.flag_s, self.flag_A,
                                                                                  self.flag_l, self.flag_C, self.flag_t,
                                                                                  self.flag_e, self.flag_F, self.flag_d)
        elif case == 'paired':
            return "/bin/bash {} {} {} {} {} {} {} {} {} {} {}".format(self.name, self.flag_s, self.flag_A,
                                                                       self.flag_l, self.flag_r, self.flag_C,
                                                                       self.flag_t, self.flag_P, self.flag_e,
                                                                       self.flag_F, self.flag_d)
        else:
            raise ValueError("unknown case")

    def __repr__(self, case: str = 'paired'):
        if case == 'single':
            return "/bin/bash {} {} {} {} -r null {} {} -P false {} {} {}".format(self.name, self.flag_s, self.flag_A,
                                                                                  self.flag_l, self.flag_C, self.flag_t,
                                                                                  self.flag_e, self.flag_F, self.flag_d)
        elif case == 'paired':
            return "/bin/bash {} {} {} {} {} {} {} {} {} {} {}".format(self.name, self.flag_s, self.flag_A,
                                                                       self.flag_l, self.flag_r, self.flag_C,
                                                                       self.flag_t, self.flag_P, self.flag_e,
                                                                       self.flag_F, self.flag_d)
        else:
            raise ValueError("unknown case")


class DeliverHaplotyperVC(Script):
    """
    Constructs the deliver_haplotyperVC.sh commands.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s output"
        self.flag_r = '-r Inputs/somaticvariants.vcf.gz'
        self.flag_j = "-j Jsons/SomaticMasterWorkflow.FilledIn.json"
        self.flag_f = '-f Delivery'  # for iforge testing
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/deliver_haplotyperVC.sh'.format(self.shell_path)
        self.type = 'deliver_haplotyperVC.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_r, self.flag_j, self.flag_f, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_r, self.flag_j, self.flag_f, self.flag_F, self.flag_d)


class BQSR(Script):
    """
    Constructs the bqsr.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_S = '-S /usr/local/apps/bioapps/sentieon/sentieon-genomics-201808'  # for iforge testing
        self.flag_G = "-G Reference/Homo_sapiens_assembly38.fasta"
        self.flag_t = '-t 40'
        self.flag_b = '-b Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_k = "-k Reference/Mills_and_1000G_gold_standard.inders.hg38.vcf"
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/bqsr.sh'.format(self.shell_path)
        self.type = 'bqsr.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_G, self.flag_t, self.flag_b, self.flag_k,
                   self.flag_e, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_G, self.flag_t, self.flag_b, self.flag_k,
                   self.flag_e, self.flag_F, self.flag_d)


class VQSR(Script):
    """
    Constructs the vqsr.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_S = '-S /usr/local/apps/bioapps/sentieon/sentieon-genomics-201808'  # for iforge testing
        self.flag_G = "-G Reference/Homo_sapiens_assembly38.fasta"
        self.flag_t = '-t 40'
        self.flag_V = '-V Inputs/somaticvariants.vcf.gz'
        self.flag_r = '-r \"\'--resource /projects/bioinformatics/DataPacks/human/gatk_bundle_Oct_2017/' \
                      'gatk_bundle_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource_param 1000G,' \
                      'known=false,training=true,truth=false,prior=10.0 --resource /projects/bioinformatics/' \
                      'DataPacks/human/gatk_bundle_Oct_2017/gatk_bundle_hg38/1000G_omni2.5.hg38.vcf.gz ' \
                      '--resource_param omni,known=false,training=true,truth=false,prior=12.0 --resource /projects/' \
                      'bioinformatics/jallen17/Reference/dbsnp_138.hg38.vcf --resource_param dbsnp,known=true,' \
                      'training=false,truth=false,prior=2.0 --resource /projects/bioinformatics/DataPacks/human/' \
                      'gatk_bundle_Oct_2017/gatk_bundle_hg38/hapmap_3.3.hg38.vcf.gz --resource_param hapmap,known=' \
                      'false,training=true,truth=true,prior=15.0\'\"'
        self.flag_R = '-R \"\'--resource /projects/bioinformatics/jallen17/Reference/dbsnp_138.hg38.vcf ' \
                      '--resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 --resource /projects/' \
                      'bioinformatics/jallen17/Reference/Mills_and_1000G_gold_standard.indels.hg38.vcf ' \
                      '--resource_param Mills,known=false,training=true,truth=true,prior=12.0\'\"'
        self.flag_a = '-a \"\'--annotation DP --annotation QD --annotation FS --annotation SOR --annotation MQ ' \
                      '--annotation MQRankSum --annotation ReadPosRankSum \'\"'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/vqsr.sh'.format(self.shell_path)
        self.type = 'vqsr.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_G, self.flag_t, self.flag_b, self.flag_k,
                   self.flag_e, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_G, self.flag_t, self.flag_b, self.flag_k,
                   self.flag_e, self.flag_F, self.flag_d)


class Alignment(Script):
    """
    Constructs the alignment.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_L = "-L fake_lib"
        self.flag_f = "-f normal"
        self.flag_c = "-c NCSA"
        self.flag_l = '-l Inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz'
        self.flag_r = '-r Inputs/WGS_chr1_5X_E0.005_L1_read2.fastq.gz'
        self.flag_G = "-G Reference/Homo_sapiens_assembly38.fasta"
        self.flag_K = "-K 10000000"
        self.flag_o = "'-M'"
        self.flag_S = '-S /usr/local/apps/bioapps/python/Python-3.6.1/bin' # for iforge testing
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/alignment.sh'.format(self.shell_path)
        self.type = 'alignment.sh'

    def __str__(self, case: str = 'paired'):
        if case == 'single':
            return "/bin/bash {} {} {} {} {} {} {} -r null {} {} {} {} {} -P false {} {} {}".\
                format(self.name, self.flag_s, self.flag_p, self.flag_L, self.flag_f, self.flag_c, self.flag_l,
                       self.flag_G, self.flag_K, self.flag_o, self.flag_S, self.flag_t,
                       self.flag_e, self.flag_F, self.flag_d)
        elif case == 'paired':
            return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".\
                format(self.name, self.flag_s, self.flag_p, self.flag_L, self.flag_f, self.flag_c, self.flag_l,
                       self.flag_r, self.flag_G, self.flag_K, self.flag_o, self.flag_S, self.flag_t, self.flag_P,
                       self.flag_e, self.flag_F, self.flag_d)
        else:
            raise ValueError("unknown case")

    def __repr__(self, case: str = 'paired'):
        if case == 'single':
            return "/bin/bash {} {} {} {} {} {} {} -r null {} {} {} {} {} -P false {} {} {}".\
                format(self.name, self.flag_s, self.flag_p, self.flag_L, self.flag_f, self.flag_c, self.flag_l,
                       self.flag_G, self.flag_K, self.flag_o, self.flag_S, self.flag_t,
                       self.flag_e, self.flag_F, self.flag_d)
        elif case == 'paired':
            return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".\
                format(self.name, self.flag_s, self.flag_p, self.flag_L, self.flag_f, self.flag_c, self.flag_l,
                       self.flag_r, self.flag_G, self.flag_K, self.flag_o, self.flag_S, self.flag_t, self.flag_P,
                       self.flag_e, self.flag_F, self.flag_d)
        else:
            raise ValueError("unknown case")


class DeDup(Script):
    """
    Constructs the dedup.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_b = '-b Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_S = '-S /usr/local/apps/bioapps/sentieon/sentieon-genomics-201808'  # for iforge testing
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/bqsr.sh'.format(self.shell_path)
        self.type = 'bqsr.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_t, self.flag_e, self.flag_S, self.flag_F,
                   self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_t, self.flag_e, self.flag_S, self.flag_F,
                   self.flag_d)


class DeliverAlignment(Script):
    """
    Constructs the deliver_alignment.sh commands.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s output"
        self.flag_b = '-b Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_j = "-j Jsons/Runalignment.FilledIn.json"
        self.flag_f = '-f Delivery'  # for iforge testing
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/deliver_alignment.sh'.format(self.shell_path)
        self.type = 'deliver_alignment.sh'

    def __str__(self, case: str = 'paired'):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_j, self.flag_f, self.flag_F, self.flag_d)

    def __repr__(self, case: str = 'paired'):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_j, self.flag_f, self.flag_F, self.flag_d)


class MergeBams(Script):
    """
    Constructs the merge_bams.sh commands.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_b = '-b Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_S = '-S /usr/local/apps/bioapps/sentieon/sentieon-genomics-201808'  # for iforge testing
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.flie'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/merge_bams.sh'.format(self.shell_path)
        self.type = 'merge_bams.sh'

    def __str__(self, case: str = 'paired'):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_j, self.flag_f, self.flag_F, self.flag_d)

    def __repr__(self, case: str = 'paired'):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_j, self.flag_f, self.flag_F, self.flag_d)


class Mutect(Script):
    """
    Constructs the mutect.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_N = '-N Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_T = '-T Inputs/WGS_chr20_21_22_tumor.bam'
        self.flag_g = "-g Reference/Homo_sapiens_assembly38.fasta"
        self.flag_G = '-G /usr/local/apps/bioapps/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef'
        self.flag_J = '-J /usr/local/apps/bioapps/java/java-1.8-64bit/bin'
        self.flag_j = '-j \"\'-Xms2G -Xmx8G\'\"'
        self.flag_B = '-B /usr/local/apps/bioapps/bcftools/bcftools-1.5'
        self.flag_Z = '-Z /usr/local/apps/bioapps/bcftools/htslib-1.3.1/bin'
        self.flag_S = '-S /usr/local/apps/bioapps/samtools/samtools-1.5'
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_D = '-D {]/../perl/fixDP.pl'.format(self.shell_path)
        self.flag_o = '-o \"\'--dbsnp /projects/bioinformatics/jallen1 /Reference/dbsnp_138.hg38.vcf\'\"'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/mutect.sh'.format(self.shell_path)
        self.type = 'mutect.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_N, self.flag_T, self.flag_g, self.flag_G, self.flag_J,
                   self.flag_j, self.flag_B, self.flag_Z, self.flag_S, self.flag_t, self.flag_e, self.flag_D,
                   self.flag_o, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_N, self.flag_T, self.flag_g, self.flag_G, self.flag_J,
                   self.flag_j, self.flag_B, self.flag_Z, self.flag_S, self.flag_t, self.flag_e, self.flag_D,
                   self.flag_o, self.flag_F, self.flag_d)


class Realignment(Script):
    """
    Constructs the realignment.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_b = '-b Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_G = "-G Reference/Homo_sapiens_assembly38.fasta"
        self.flag_k = "-k Reference/Mills_and_1000G_gold_standard.inders.hg38.vcf"
        self.flag_S = '-S /usr/local/apps/bioapps/sentieon/sentieon-genomics-201808'
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/realignment.sh'.format(self.shell_path)
        self.type = 'realignment.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_G, self.flag_k, self.flag_S, self.flag_t,
                   self.flag_e, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_b, self.flag_G, self.flag_k, self.flag_S, self.flag_t,
                   self.flag_e, self.flag_F, self.flag_d)


class Strelka(Script):
    """
    Constructs the strelka.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_N = '-N Inputs/WGS_chr20_21_22_normal.bam'
        self.flag_T = '-T Inputs/WGS_chr20_21_22_tumor.bam'
        self.flag_g = "-g Reference/Homo_sapiens_assembly38.fasta"
        self.flag_B = '-B /usr/local/apps/bioapps/bcftools/bcftools-1.5'
        self.flag_I = '-I /usr/local/apps/bioapps/strelka/strelka-2.9.2.centos6_x86_64/bin'
        self.flag_S = '-S /usr/local/apps/bioapps/samtools/samtools-1.5'
        self.flag_Z = '-Z /usr/local/apps/bioapps/bcftools/htslib-1.3.1/bin'
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_i = '-i MayomicsVC/src/perl/fixStrelka_GT_indels.pl'
        self.flag_p = '-p MayomicsVC/src/perl/fixStrelka_GT_snvs.pl'
        self.flag_o = '-o \"\'--outputCallableRegions\'\"'
        self.flag_O = '-O \"\'-m any --force-sample\'\"'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/strelka.sh'.format(self.shell_path)
        self.type = 'strelka.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_N, self.flag_T, self.flag_g, self.flag_B, self.flag_I,
                   self.flag_S, self.flag_Z, self.flag_t, self.flag_e, self.flag_i, self.flag_p, self.flag_o,
                   self.flag_O, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_N, self.flag_T, self.flag_g, self.flag_B, self.flag_I,
                   self.flag_S, self.flag_Z, self.flag_t, self.flag_e, self.flag_i, self.flag_p, self.flag_o,
                   self.flag_O, self.flag_F, self.flag_d)


class CombineVariants(Script):
    """
    Constructs the strelka.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_S = '-S Inputs/strelka.vcf.bgz'
        self.flag_T = '-T Inputs/mutect.vcf.bgz'
        self.flag_g = "-g Reference/Homo_sapiens_assembly38.fasta"
        self.flag_G = '-G /usr/local/apps/bioapps/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef'
        self.flag_J = '-J /usr/local/apps/bioapps/java/java-1.8-64bit/bin'
        self.flag_B = '-B /usr/local/apps/bioapps/bcftools/bcftools-1.5'
        self.flag_Z = '-Z /usr/local/apps/bioapps/bcftools/htslib-1.3.1/bin'
        self.flag_t = '-t 40'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_o = '-o \"\' \'\"'
        self.flag_p = '-p \"\'strelka,mutect\'\"'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/strelka.sh'.format(self.shell_path)
        self.type = 'strelka.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_T, self.flag_g, self.flag_G, self.flag_J,
                   self.flag_B, self.flag_Z, self.flag_t, self.flag_e, self.flag_o, self.flag_p, self.flag_F,
                   self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_T, self.flag_g, self.flag_G, self.flag_J,
                   self.flag_B, self.flag_Z, self.flag_t, self.flag_e, self.flag_o, self.flag_p, self.flag_F,
                   self.flag_d)


class DeliverSomaticVC(Script):
    """
    Constructs the deliver_haplotyperVC.sh commands.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s output"
        self.flag_r = '-r Inputs/somaticvariants.vcf.gz'
        self.flag_j = "-j Jsons/SomaticMasterWorkflow.FilledIn.json"
        self.flag_f = '-f Delivery'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/deliver_somaticVC.sh'.format(self.shell_path)
        self.type = 'deliver_somaticVC.sh'

    def __str__(self, case: str = 'paired'):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_r, self.flag_j, self.flag_f, self.flag_F, self.flag_d)

    def __repr__(self, case: str = 'paired'):
        return "/bin/bash {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_r, self.flag_j, self.flag_f, self.flag_F, self.flag_d)


class Haplotyper(Script):
    """
    Constructs the haplotyper.sh commands for single or paired end reads.
    """

    def __init__(self):
        Script.__init__(self)
        self.flag_s = "-s outputs/output"
        self.flag_S = '-S /usr/local/apps/bioapps/sentieon/sentieon-genomics-201808'
        self.flag_G = "-G Reference/Homo_sapiens_assembly38.fasta"
        self.flag_t = '-t 40'
        self.flag_b = '-b Inputs/sample.bam'
        self.flag_D = '-D Reference/dbsnp_138.hg38.vcf'
        self.flag_r = '-r Inputs/bqsr.recal_data.table'
        self.flag_o = '-o \"\' \'\"'
        self.flag_e = '-e Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.shell_path)
        self.flag_d = '-d'
        self.name = '{}/haplotyper.sh'.format(self.shell_path)
        self.type = 'haplotyper.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_G, self.flag_t, self.flag_b, self.flag_D,
                   self.flag_r, self.flag_o, self.flag_e, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {} {}". \
            format(self.name, self.flag_s, self.flag_S, self.flag_G, self.flag_t, self.flag_b, self.flag_D,
                   self.flag_r, self.flag_o, self.flag_e, self.flag_F, self.flag_d)


class ParameterizedTestCase(unittest.TestCase):
    """
    Test cases with parameters will inherit from this class
    Code borrowed from and adapted: https://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases
    """
    def __init__(self, methodName='runTest', param=None):
        super(ParameterizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parameterize(testcase_klass, param=None):
        """
        Create a suite containing all tests taken from the given subclass, passing them
        the parameter 'param'
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, param=param))
        return suite


class TestArgs(ParameterizedTestCase):

    def setUp(self):
        pass

    def tearDown(self):
        files = glob.glob('outputs/*')
        for f in files:
            os.remove(f)
        files = glob.glob('WGS*')
        for f in files:
            os.remove(f)

    # @unittest.skip("Testing")
    def test_no_arg(self):
        """
        Tests the script produces help output when no argument is passed
        """
        os.system("/bin/bash " + self.param.name + ' > outputs/outfile.txt')
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue('command line input: \n' in output)
        self.assertTrue("No arguments passed." in output)

    # @unittest.skip("Testing")
    def test_help_function(self):
        """
        While this theoretically works for all, the tricky part is that each script has a unique help output,
        since each takes a different set of inputs. It'll take a bit of modification to get it to compare
        the correct files. It'll also require building and maintaining a desired help-file database
        """
        os.system("/bin/bash " + self.param.name + ' -h > outputs/outfile.txt')
        desired_help = self.parse_output('{}/Verification_files/desired_help_output.txt'.format(self.param.test_path))
        output = self.parse_output('outputs/outfile.txt')
        for i in range(4, len(output)-1):
            self.assertTrue(desired_help[i-4] == output[i])

    # @unittest.skip("Testing")
    def test_nonexistent_option(self):
        """
        Test a flag that doesn't exist with a garbage test option. This should work for all as is.
        """
        os.system("/bin/bash " + self.param.name + " -Q garbage > outputs/outfile.txt")
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue('command line input: -Q garbage' in output)
        self.assertTrue("Invalid option: -Q" in output)

    # @unittest.skip("So slow")
    def test_successful_paired_end_read_and_permissions(self):
        """
        This is simply a successful run of the tool and should be generalizable. I'm including a
        check for file permissions to avoid having to duplicate code later.
        """
        os.system(self.param.__str__('paired') + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        # check that it started and ended properly
        self.assertTrue('START' in output)
        self.assertTrue("Finished trimming adapter sequences." in output)
        cutadapt_output = 'WGS_chr1_5X_E0.005_L1_read1.fastq.gz'

        # Check that it created a non-zero output
        self.assertTrue(os.path.exists(cutadapt_output) and os.path.getsize(cutadapt_output) > 0)

        # check that permissions have been set properly
        perm_check_log = subprocess.Popen(['ls', '-l', 'outputs/output.trimming.TBD.log'], stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT)
        stdout_log, stderr_log = perm_check_log.communicate()
        self.assertTrue("-rw-r-----" in str(stdout_log))
        perm_check_read1 = subprocess.Popen(['ls', '-l', 'WGS_chr1_5X_E0.005_L1_read1.fastq.gz'],
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_read1, stderr_read1 = perm_check_read1.communicate()
        self.assertTrue("-rw-r-----" in str(stdout_read1))
        perm_check_read2 = subprocess.Popen(['ls', '-l', 'WGS_chr1_5X_E0.005_L1_read2.fastq.gz'],
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_read2, stderr_read2 = perm_check_read2.communicate()
        self.assertTrue("-rw-r-----" in str(stdout_read2))

        # Test minimal permissions on the output files
        os.chmod('WGS_chr1_5X_E0.005_L1_read1.fastq.gz', 0o200)
        os.system(self.param.__str__('paired') + " > outputs/outfile.txt 2>&1 ")
        self.assertTrue(oct(os.stat('WGS_chr1_5X_E0.005_L1_read2.fastq.gz').st_mode)[-3:][1] == '4')

        os.chmod('WGS_chr1_5X_E0.005_L1_read2.fastq.gz', 0o200)
        os.system(self.param.__str__('paired') + " > outputs/outfile.txt 2>&1 ")
        self.assertTrue(oct(os.stat('WGS_chr1_5X_E0.005_L1_read2.fastq.gz').st_mode)[-3:][1] == '4')


    # @unittest.skip("So slow")
    def test_successful_single_end_read(self):
        """
        This is simply a successful run of the tool and should be generalizable
        """
        os.system(self.param.__str__('single') + " > outputs/outfile.txt 2>&1")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue('START' in output)
        self.assertTrue("Finished trimming adapter sequences." in output)
        cutadapt_log = 'outputs/output.cutadapt.log'
        self.assertTrue(os.path.exists(cutadapt_log) and os.path.getsize(cutadapt_log) > 0)

    # @unittest.skip("So slow")
    def test_read_flags_with_bad_input(self):
        """
        Most of the scripts have some sort of input, so this will probably be generalizable to a degree.
        It simply tries some dummy/garbage files for read inputs to make sure the tool isn't trying to align
        text that isn't genomic.
        """
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")

        # test left and right read flags
        flags_to_test = ['flag_l', 'flag_r']
        garbage_test_files = {'dummy_test_blank.fastq':
                              "garbage_test_files/dummy_test_blank.fastq is empty or does not exist.",
                              'dummy_test_text.fastq':
                                  "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found "
                                  "'Lorem ipsu'",
                              'dummy_test_text_with_at.fastq':
                                  "cutadapt: error: Line 3 in FASTQ file is expected to "
                                  "start with '+', but found 'Suspendiss'",
                              'WGS_chr1_5X_E0.005_L1_read1.fastq.':
                                  'WGS_chr1_5X_E0.005_L1_read1.fastq. is empty or does not exist'}
        for flag in flags_to_test:
            for garbage_test in garbage_test_files.keys():
                temp_flag = copy.deepcopy(self.param.__dict__[flag])
                manip_flag = self.param.__dict__[flag]
                if "dummy" in garbage_test:
                    manip_flag = manip_flag.split(' ')[0] + ' {}/garbage_test_files/'.format(self.param.test_path) \
                                 + garbage_test
                else:
                    manip_flag = manip_flag.split(' ')[0] + ' Inputs/' + garbage_test
                self.param.__dict__[flag] = manip_flag
                os.system(str(self.param) + " > outputs/outfile.txt 2>&1 ")
                output = self.parse_output('outputs/output.trimming.TBD.log')
                log = self.parse_output('outputs/output.cutadapt.log')
                output = ''.join(output)
                log = ''.join(log)
                if 'Cutadapt Read 1 and 2 failure' in output:
                    self.assertTrue(garbage_test_files[garbage_test] in log)
                else:
                    self.assertTrue(garbage_test_files[garbage_test] in output)
                self.param.__dict__[flag] = temp_flag
                try:
                    os.remove(garbage_test)
                except OSError:
                    pass

    # @unittest.skip("So slow")
    def test_garbage_adapters(self):
        """
        This tests trim_sequences call for the adapter files. This may be generalizabale, since other scripts
        will call in outside files as well.
        """
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")

        tests = {'dummy_test_blank.fastq': "garbage_test_files/dummy_test_blank.fastq is empty or does not exist.",
                 'dummy_test_text.fastq': "At line 1: Expected '>' at beginning of FASTA record, but got 'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'",
                 'dummy_test_text_with_gt.fastq': "is not a valid IUPAC code. Use only characters XACGTURYSWKMBDHVN.",
                 'TruSeqAdapters.fasta': "TruSeqAdapters.fasta is empty or does not exist"}
        for test in tests.keys():
            temp_flag = copy.deepcopy(self.param.__dict__['flag_A'])
            manip_flag = self.param.__dict__['flag_A']
            manip_flag = manip_flag.split(' ')[0] + ' {}/garbage_test_files/'.format(self.param.test_path) + test
            self.param.__dict__['flag_A'] = manip_flag
            os.system(self.param.__str__('paired') + " > outputs/outfile.txt 2>&1 ")
            output = self.parse_output('outputs/output.trimming.TBD.log')
            log = self.parse_output('outputs/output.cutadapt.log')
            output = ''.join(output)
            log = ''.join(log)
            with self.subTest(test=test):
                if 'Cutadapt Read 1 and 2 failure' in output:
                    self.assertTrue(tests[test] in log)
                else:
                    self.assertTrue(tests[test] in output)
            self.param.__dict__['flag_A'] = temp_flag
            try:
                os.remove(test)
            except OSError:
                pass

    @unittest.skip("So slow")
    def test_bad_env_file(self):
        """
        This simply uses a non-existant environmental file to test that the script checks for this.
        """
        tests = {'envprof_fake.file': "No such file or directory"}
        for test in tests.keys():
            temp_flag = copy.deepcopy(self.param.__dict__['flag_e'])
            manip_flag = self.param.__dict__['flag_e']
            manip_flag = manip_flag.split(' ')[0] + ' ' + test
            self.param.__dict__['flag_e'] = manip_flag
            os.system(self.param.__str__('paired') + " > outputs/outfile.txt 2>&1 ")
            output = self.parse_output('outputs/outfile.txt')
            output = ''.join(output)
            self.assertTrue(tests[test] in output)
            self.param.__dict__['flag_e'] = temp_flag

    @unittest.skip("So slow")
    def test_bad_cutadapt_path(self):
        """
        Trim_sequences is the only one that uses cutadapt, though I may be able to generalize this script
        to 'test bad tool path' once I get rolling on the others.
        """
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")

        # Test bad cutadapt path
        os.system("/bin/bash {} {} {} {} {} -C /usr/fake {} {} {} {} {} > outputs/outfile.txt 2>&1".
                  format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_l, self.param.flag_r,
                         self.param.flag_t, self.param.flag_P, self.param.flag_e, self.param.flag_F, self.param.flag_d))
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue("REASON=Cutadapt directory /usr/fake is not a directory or does not exist." in output)

    @unittest.skip("So slow")
    def test_bad_thread_options(self):
        """
        This tests trim_sequences thread option. It should return errors for having too high a thread count,
        though this may vary by machine.
        """
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")

        values = [321, 3299, 12322]
        for number in values:
            os.system("/bin/bash {} {} {} {} {} {} {} {} {} {}".\
                format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_l, self.param.flag_r,
                       self.param.flag_C, self.param.flag_P, self.param.flag_e, self.param.flag_F, self.param.flag_d)
                      + " -t " + str(number) + " > outputs/outfile.txt 2>&1")
            output = self.parse_output('outputs/output.trimming.TBD.log')
            log = self.parse_output('outputs/output.cutadapt.log')
            output = ''.join(output)
            log = ''.join(log)
            if 'Finished trimming adapter sequences.' in output:
                self.assertTrue(True)
            else:
                self.assertTrue('Cutadapt Read 1 and 2 failure.' in output)

    @unittest.skip("Testing")
    def test_paired_options(self):
        """
        Not every script will have a true/false/read value to test; this only works with ones that do
        """
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")

        tests = ['True', 'T', 'False', "F"]
        for test in tests:
            os.system("/bin/bash {} {} {} {} {} {} {} {} {} {}".\
                format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_l, self.param.flag_r,
                       self.param.flag_C, self.param.flag_t, self.param.flag_e, self.param.flag_F, self.param.flag_d) +
                      ' -P ' + test + " > outputs/outfile.txt 2>&1 ")
            output = self.parse_output('outputs/output.trimming.TBD.log')
            output = ''.join(output)
            self.assertTrue("REASON=Incorrect argument for paired-end option -P. Must be set to true or false."
                            in output)

    @unittest.skip("Testing")
    def test_incorrect_read_options(self):
        """
        Not every script will have a true/false/read value to test; this only works with ones that do
        """
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")
        os.system("/bin/bash {} {} {} {} {} {} {} -P false {} {} {}". \
                  format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_l, self.param.flag_r,
                         self.param.flag_C, self.param.flag_t, self.param.flag_e,
                         self.param.flag_F, self.param.flag_d) + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue("REASON=User specified Single End option, but did not set read 2 option -r to null." in output)
        os.system("/bin/bash {} {} {} {} {} {} {} -P true {} {}". \
                  format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_l,
                         self.param.flag_C, self.param.flag_t, self.param.flag_e,
                         self.param.flag_F, self.param.flag_d) + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue("REASON=Missing read 2 option: -r. If running a single-end job, set -r null in command." in output)
        os.system("/bin/bash {} {} {} -l null {} {} {} -P false {} {} {}". \
                  format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_r,
                         self.param.flag_C, self.param.flag_t, self.param.flag_e,
                         self.param.flag_F, self.param.flag_d) + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue("REASON=Input read 1 file null is empty or does not exist." in output)
        os.system("/bin/bash {} {} {} -r null {} {} {} {} -P true {} {}". \
                  format(self.param.name, self.param.flag_s, self.param.flag_A, self.param.flag_l,
                         self.param.flag_C, self.param.flag_t, self.param.flag_e,
                         self.param.flag_F, self.param.flag_d) + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue("REASON=Input read 2 file null is empty or does not exist." in output)

    @unittest.skip('Test')
    def test_missing_option_values(self):
        """
        Should work with any of the scripts. Note that -d flag is ommitted since it does not have a value
        already
        """
        attributes = list(self.param.__dict__.keys())
        attributes.remove('flag_d')
        options = list([a for a in attributes if "flag" in a])
        for flag in options:
            temp_flag = copy.deepcopy(self.param.__dict__[flag])
            manip_flag = self.param.__dict__[flag]
            manip_flag = manip_flag.split(' ')[0]
            self.param.__dict__[flag] = manip_flag
            os.system(str(self.param) + " > outputs/outfile.txt 2>&1 ")
            output = self.parse_output('outputs/outfile.txt')
            output = ''.join(output)
            self.assertTrue("Error with option " + manip_flag + " in command. Option passed incorrectly or without argument." in output)
            self.param.__dict__[flag] = temp_flag

    @unittest.skip('testing')
    def test_file_permissions(self):
        """
        Should work with any of the scripts
        """
        os.chmod('Inputs', 0o000)
        os.system(str(self.param) + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue('is empty or does not exist' in output)
        os.chmod('Inputs', 0o755)

        os.system('chmod 000 outputs')
        os.system(str(self.param) + " > outfile.txt 2>&1 ")
        output = self.parse_output('outfile.txt')
        output = ''.join(output)
        self.assertTrue('Permission denied' in output)
        os.remove('outfile.txt')
        os.chmod('outputs', 0o755)

        os.chmod(self.param.shell_path, 0o000)
        os.system(str(self.param) + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue('Permission denied' in output)
        os.chmod(self.param.shell_path, 0o755)

    @unittest.skip("Testing")
    def test_logs_are_truncated(self):
        # first run creates logs
        os.system("/bin/bash {} -s outputs/output -A garbage_test_files/dummy_test_text.fastq {} {} {} -P true {} {} {}"
                  " {}".format(self.param.name, self.param.flag_r, self.param.flag_l, self.param.flag_C,
                               self.param.flag_t, self.param.flag_e, self.param.flag_F, self.param.flag_d) +
                  " > outputs/outfile.txt 2>&1 ")
        output_stdout = self.parse_output('outputs/output.trimming.TBD.log')
        output_stdout_test = output_stdout[-2:]
        output_stdout_test = ''.join(output_stdout_test)
        output_stdout = ''.join(output_stdout)
        output_cutlog = self.parse_output('outputs/output.cutadapt.log')
        output_cutlog_test = output_cutlog[-2:]
        output_cutlog_test = ''.join(output_cutlog_test)
        output_cutlog = ''.join(output_cutlog)

        # second run
        time.sleep(2)
        os.system("/bin/bash {} -s outputs/output -A garbage_test_files/dummy_test_text_with_gt.fastq {} {} {} -P true"
                  " {} {} {} {}".format(self.param.name, self.param.flag_r, self.param.flag_l, self.param.flag_C,
                                        self.param.flag_t, self.param.flag_e, self.param.flag_F, self.param.flag_d) +
                  " > outputs/outfile.txt 2>&1 ")
        output_stdout2 = self.parse_output('outputs/output.trimming.TBD.log')
        output_stdout2 = ''.join(output_stdout2)
        output_cutlog2 = self.parse_output('outputs/output.cutadapt.log')
        output_cutlog2 = ''.join(output_cutlog2)

        # The logs should be different and the second log shouldn't be contained in the first
        self.assertNotEqual(output_stdout, output_stdout2)
        self.assertNotEqual(output_cutlog, output_cutlog2)
        self.assertTrue(output_stdout_test not in output_stdout2)
        self.assertTrue(output_cutlog_test not in output_cutlog2)

    @staticmethod
    def parse_output(file):
        output = []
        for line in open(file, 'r'):
            output.append(line)
        return output


if __name__ == "__main__":
    scripts = ["trim_sequences.sh", 'deliver_haplotyperVC.sh', 'bqsr.sh', 'vqsr.sh', 'alignment.sh', 'dedup.sh',
               'deliver_alignment.sh', 'merge_bams.sh', 'mutect.sh', 'realignment.sh', 'strelka.sh',
               'combine_variants.sh', 'deliver_somaticVC.sh', 'haplotyper.sh']
    try:
        idx = scripts.index(sys.argv[1])
    except ValueError:
        print("Argument must be the script to test and the output_file/log_name to use.")
    if idx == 0:
        test_script = Trimming()
    elif idx == 1:
        test_script = DeliverHaplotyperVC()
    elif idx == 2:
        test_script = BQSR()
    elif idx == 3:
        test_script = VQSR()
    elif idx == 4:
        test_script = Alignment()
    elif idx == 5:
        test_script = DeDup()
    elif idx == 6:
        test_script = DeliverAlignment()
    elif idx == 7:
        test_script = MergeBams()
    elif idx == 8:
        test_script = Mutect()
    elif idx == 9:
        test_script = Realignment()
    elif idx == 10:
        test_script = Strelka()
    elif idx == 11:
        test_script = CombineVariants()
    elif idx == 12:
        test_script = DeliverSomaticVC()
    elif idx == 13:
        test_script = Haplotyper()

    suite = unittest.TestSuite()
    suite.addTest(ParameterizedTestCase.parameterize(TestArgs, param=test_script))

    unittest.TextTestRunner(verbosity=2).run(suite)



