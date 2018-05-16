import ToolPath_Parser
import unittest


class TestParsingTools(unittest.TestCase):

    testInputLines = ['#List of Tools', 'BWA="/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']
    testDummyInputLine = ['#List of Tools', 'BWA="/usr/local/apps====/bioapps=======/bwa/bwa-0.7.16/bwa"']
    
    testList = ['BWA="/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']  
    testEqualSign = ['BWA"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']
   
    testListofList = [['BWA', '"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']] 
    testEmptyExec = [['BWA', '']]
    testDoubleQuotes = [['BWA', '/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa']] 
    testSpecialChar = [['BWA', '"/us!r/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']]
    testToolRepeat = [['BWA','"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"'], ['BWA', '"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']]

    testTools = ['BWA']
    testPaths = ['"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']
    testDict = {'"CallReadMappingTask.ReadMappingTask.BWA"': '"File"'}


    def testCommentRemoval(self):
        self.assertEqual(ToolPath_Parser.commentRemoval(self.testInputLines), ['BWA="/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"'])
 

    def testEqualSignAbsence(self):
        with self.assertRaises(SystemExit):
            ToolPath_Parser.keyValuePairCreation(self.testEqualSign)
    

    def testKeyValuePairCreation(self):
        self.assertEqual(ToolPath_Parser.keyValuePairCreation(self.testList), [['BWA', '"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']] )


    def testExecutableAbsense(self):
        with self.assertRaises(SystemExit):
            ToolPath_Parser.toolsandPaths(self.testEmptyExec)


    def testDoubleQuotesAbsense(self):
        with self.assertRaises(SystemExit):
            ToolPath_Parser.toolsandPaths(self.testDoubleQuotes)


    def testSpecialCharPresence(self):
        with self.assertRaises(SystemExit):
            ToolPath_Parser.toolsandPaths(self.testSpecialChar)


    def testToolRepetition(self):
        with self.assertRaises(SystemExit):
            ToolPath_Parser.toolsandPaths(self.testToolRepeat)


    def test_ToolsandPaths(self):
        self.assertEqual(ToolPath_Parser.toolsandPaths(self.testListofList), (['BWA'], ['"/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa"']))


    def test_ExecutablesCapture(self):
        self.assertEqual(ToolPath_Parser.executablesCapture(self.testDict, self.testTools, self.testPaths), {'"CallReadMappingTask.ReadMappingTask.BWA"': '/usr/local/apps/bioapps/bwa/bwa-0.7.16/bwa'})


if __name__ == '__main__':
    unittest.main()
