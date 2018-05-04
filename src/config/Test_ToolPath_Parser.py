import ToolPath_Parser
import unittest


class TestParsingTools(unittest.TestCase):

    testInputLines = ['#List of Tools', 'BWA="/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']
    testDummyInputLine = ['#List of Tools', 'BWA="/usr/local/apps====/bioapps=======/bwa/bwa-0.7.16/bwa"']
    testList = ['BWA="/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']  
    testListofList = [['BWA', '"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']]
    testTools = ['BWA']
    testPaths = ['"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']
    testDict = {'"CallReadMappingTask.ReadMappingTask.BWA"': '"File"'}


    def testCommentRemoval(self):
        self.assertEqual(ToolPath_Parser.commentRemoval(self.testInputLines), ['BWA="/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"'])
  

    def testKeyValuePairCreation(self):
        self.assertEqual(ToolPath_Parser.keyValuePairCreation(self.testList), [['BWA', '"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']] )


    def test_ToolsandPaths(self):
        self.assertEqual(ToolPath_Parser.toolsandPaths(self.testListofList), (['BWA'], ['"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']))


    def test_ExecutablesCapture(self):
        self.assertEqual(ToolPath_Parser.executablesCapture(self.testDict, self.testTools, self.testPaths), {'"CallReadMappingTask.ReadMappingTask.BWA"': '/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa'})


if __name__ == '__main__':
    unittest.main()
