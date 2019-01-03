import sys
import os
import subprocess
import unittest
from subprocess import Popen, PIPE
import glob
import copy


class Script:
    """
    Trim_sequences is done, and the other scripts will have similar setups that have the flags as attributes so tests
    can step through the list of attributes.
    """
    __cwd = os.getcwd()

    def __init__(self):
        self.path = '{}/../../src/shell'.format(self.__cwd)


class Trimming(Script):
    """
    Runs trim_sequences.sh. Contains the parts to run the code. Each attribute represents a particular flag, that way
    we can step through the flags and perform tests on each.

    TODO: os.path to make tests more generic
    """

    def __init__(self, output, threads):
        Script.__init__(self)
        self.flag_s = "-s outputs/{}".format(output)
        self.flag_A = '-A ../../../Inputs/TruSeqAdaptors.fasta'
        self.flag_l = '-l ../../../Inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz'
        self.flag_r = '-r ../../../Inputs/WGS_chr1_5X_E0.005_L1_read2.fastq.gz'
        # self.flag_C = '-C /usr/local/apps/bioapps/python/Python-3.6.1/bin' # for iforge testing
        self.flag_C = '-C /usr/bin' # for local testing
        self.flag_t = '-t {}'.format(threads)
        self.flag_P = '-P true'
        self.flag_e = '-e ../../../Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.path)
        self.flag_d = '-d'
        self.name = '{}/trim_sequences.sh'.format(self.path)
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

class Alignment(Script):

    def __init__(self, output, threads):
        Script.__init__(self)
        self.flag_s = "-s outputs/{}".format(output)
        self.flag_A = '-A ../../../Inputs/TruSeqAdaptors.fasta'
        self.flag_l = '-l ../../../Inputs/WGS_chr1_5X_E0.005_L1_read1.fastq.gz'
        self.flag_r = '-r ../../../Inputs/WGS_chr1_5X_E0.005_L1_read2.fastq.gz'
        # self.flag_C = '-C /usr/local/apps/bioapps/python/Python-3.6.1/bin' # for iforge testing
        self.flag_C = '-C /usr/bin' # for local testing
        self.flag_t = '-t {}'.format(threads)
        self.flag_P = '-P true'
        self.flag_e = '-e ../../../Config/EnvProfile.file'
        self.flag_F = '-F {}/shared_functions.sh'.format(self.path)
        self.flag_d = '-d'
        self.name = '{}/trim_sequences.sh'.format(self.path)
        self.type = 'alignment.sh'

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


class MergeBams(Script):
    pass


class DeDup(Script):
    pass


class Realignment(Script):
    pass


class BQSR(Script):
    pass


class Haplotyper(Script):
    pass


class VQSR(Script):
    pass


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

    def test_cutadapt_installed(self):
        os.system("apt list cutadapt > outputs/outfile.txt 2>&1")
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue("[installed]" in output)

    def test_no_arg(self):
        os.system("/bin/bash " + self.param.name + ' > outputs/outfile.txt')
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue('command line input: \n' in output)
        self.assertTrue("No arguments passed." in output)

    def test_help_function(self):
        os.system("/bin/bash " + self.param.name + ' -h > outputs/outfile.txt')
        desired_help = self.parse_output('Verification_files/desired_help_output.txt')
        output = self.parse_output('outputs/outfile.txt')
        for i in range(4, len(output)-1):
            self.assertTrue(desired_help[i-4] == output[i])

    def test_nonexistent_option(self):
        os.system("/bin/bash " + self.param.name + " -Q garbage > outputs/outfile.txt")
        output = self.parse_output('outputs/outfile.txt')
        output = ''.join(output)
        self.assertTrue('command line input: -Q garbage' in output)
        self.assertTrue("Invalid option: -Q" in output)

    @unittest.skip("So slow")
    def test_successful_paired_end_read(self):
        os.system(self.param.__str__('paired') + " > outputs/outfile.txt 2>&1 ")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue('START' in output)
        self.assertTrue("Finished trimming adapter sequences." in output)
        cutadapt_log = 'outputs/output.cutadapt.log'
        self.assertTrue(os.path.exists(cutadapt_log) and os.path.getsize(cutadapt_log) > 0)

    @unittest.skip("So slow")
    def test_successful_single_end_read(self):
        os.system(self.param.__str__('single') + " > outputs/outfile.txt 2>&1")
        output = self.parse_output('outputs/output.trimming.TBD.log')
        output = ''.join(output)
        self.assertTrue('START' in output)
        self.assertTrue("Finished trimming adapter sequences." in output)
        cutadapt_log = 'outputs/output.cutadapt.log'
        self.assertTrue(os.path.exists(cutadapt_log) and os.path.getsize(cutadapt_log) > 0)

    def test_options_for_trimming(self):
        if self.param.type != 'trim_sequences.sh':
            print("Only valid for trim sequences")
            return unittest.skip("Only valid for trim_sequences")
        options = ['flag_l', 'flag_r']
        tests = {'dummy_test_blank.fastq': "file garbage_test_files/dummy_test_blank.fastq is empty or does not exist.",
                 'dummy_test_text.fastq': "cutadapt: error: Line 1 in FASTQ file is expected to start with '@', but found 'Lorem ipsu'",
                 'dummy_test_text_with_at.fastq': "cutadapt: error: Line 3 in FASTQ file is expected to start with '+', but found 'Suspendiss'"}
        for option in options:
            for test in tests.keys():
                temp_flag = copy.deepcopy(self.param.__dict__[option])
                manip_flag = self.param.__dict__[option]
                manip_flag = manip_flag.split(' ')[0] + ' garbage_test_files/' + test
                self.param.__dict__[option] = manip_flag
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
                self.param.__dict__[option] = temp_flag

        options = ['flag_A']
        tests = {'dummy_test_blank.fastq': "file garbage_test_files/dummy_test_blank.fastq is empty or does not exist.",
                 'dummy_test_text.fastq': "At line 1: Expected '>' at beginning of FASTA record, but got 'Lorem ipsum dolor sit amet, consectetur adipiscing elit.'",
                 'dummy_test_text_with_gt.fastq': "is not a valid IUPAC code. Use only characters XACGTURYSWKMBDHVN."}
        for option in options:
            for test in tests.keys():
                temp_flag = copy.deepcopy(self.param.__dict__[option])
                manip_flag = self.param.__dict__[option]
                manip_flag = manip_flag.split(' ')[0] + ' garbage_test_files/' + test
                self.param.__dict__[option] = manip_flag
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
                self.param.__dict__[option] = temp_flag

    @staticmethod
    def parse_output(file):
        output = []
        for line in open(file, 'r'):
            output.append(line)
        return output


if __name__ == "__main__":
    scripts = ["trim_sequences.sh", 'alignment.sh', 'merge_bams.sh', 'dedup.sh',
               'realignment.sh', 'bqsr.sh', 'haplotyper.sh', 'vqsr.sh']
    try:
        idx = scripts.index(sys.argv[1])
    except ValueError:
        print("Argument must be the script to test and the output_file/log_name to use.")
    if idx == 0:
        test_script = Trimming('output', 0)
    elif idx == 1:
        test_script = Alignment('output', 20)
    elif idx == 2:
        test_script = MergeBams()
    elif idx == 3:
        test_script = DeDup()
    elif idx == 4:
        test_script = Realignment()
    elif idx == 5:
        test_script = BQSR()
    elif idx == 6:
        test_script = Haplotyper()
    elif idx == 7:
        test_script = VQSR()


    suite = unittest.TestSuite()
    suite.addTest(ParameterizedTestCase.parameterize(TestArgs, param=test_script))

    unittest.TextTestRunner(verbosity=2).run(suite)
