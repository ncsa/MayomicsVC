import sys
import os
import subprocess
import unittest
from subprocess import Popen, PIPE


class Script:
    __cwd = os.getcwd()

    def __init__(self, flags: list):
        self.flags = flags
        self.path = '{}/../../src/shell'.format(self.__cwd)

    def __str__(self):
        return "/bin/bash {}".format(self.__class__)

    def __repr__(self):
        return "/bin/bash {}".format(self.__class__)


class Trimming(Script):
    __cwd = os.getcwd()

    def __init__(self, threads, paired):
        flags = ["s", "A", "l", "r", "C", "t", "P", "e", "F", "d"]
        Script.__init__(self, flags)
        self.flag_s = "-s {}/../../../Output/testrun".format(self.__cwd)
        self.flag_A = '-A {}/../../../Inputs/TruSeqAdaptors.fasta'.format(self.__cwd)
        self.flag_l = '-l {}/../../../Inputs/WGS_chr1_5X_E0.005_L1_read1.fastq'.format(self.__cwd)
        self.flag_r = '-r {}/../../../Inputs/WGS_chr1_5X_E0.005_L1_read2.fastq'.format(self.__cwd)
        self.flag_C = '-C /usr/local/apps/bioapps/python/Python-3.6.1/bin'
        self.flag_t = '-t {}'.format(threads)
        self.flag_P = '-P {}'.format(paired)
        self.flag_e = '-e {}/../../../Inputs/env.file'.format(self.__cwd)
        self.flag_F = '-F {}/shared_functions.sh'.format(self.path)
        self.flag_d = ''
        self.name = '../../src/shell/trim_sequences.sh'

    def __str__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {}".format(self.name, self.flag_s, self.flag_A, self.flag_r,
                                                                   self.flag_l, self.flag_C, self.flag_t, self.flag_P,
                                                                   self.flag_e, self.flag_F, self.flag_d)

    def __repr__(self):
        return "/bin/bash {} {} {} {} {} {} {} {} {} {} {}".format(self.name, self.flag_s, self.flag_A, self.flag_r,
                                                                   self.flag_l, self.flag_C, self.flag_t, self.flag_P,
                                                                   self.flag_e, self.flag_F, self.flag_d)

class Alignment(Script):
    __cwd = os.getcwd()

    pass


class MergeBams(Script):
    __cwd = os.getcwd()
    pass


class DeDup(Script):
    __cwd = os.getcwd()
    pass


class Realignment(Script):
    __cwd = os.getcwd()
    pass


class BQSR(Script):
    __cwd = os.getcwd()
    pass


class Haplotyper(Script):
    __cwd = os.getcwd()
    pass


class VQSR(Script):
    __cwd = os.getcwd()
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
        # os.remove('outfile.txt')
        pass

    def test_no_arg(self):
        os.system("/bin/bash " + self.param.name + '> outfile.txt')
        output = self.parse_output('outfile.txt')
        self.assertTrue('command line input: \n' in output[4])
        self.assertTrue("No arguments passed." in output[7])

    def test_help_function(self):
        os.system("/bin/bash " + self.param.name + ' -h > outfile.txt')
        output = self.parse_output('outfile.txt')
        self.assertTrue('command line input: -h' in output[4])

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
        idx = scripts.index(sys.argv[1]+".sh")
    except ValueError:
        print("Argument must be the script to test and the output_file/log_name to use.")
    if idx == 0:
        test_script = Trimming(20, 'true')
    elif idx == 1:
        test_script = Alignment()
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
