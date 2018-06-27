import src.config.parsing.parser as parse
import unittest


class TestParsingTools(unittest.TestCase):

    # Create an instance of the Parser class
    parser = parse.Parser(job_id="NA")

    def test_remove_comments(self):
        input_lines = ["# Comment line", "      # Whitespace with comment", 'Key="Value"']
        filtered_lines = parse.Parser.remove_comments(input_lines)
        self.assertEqual(filtered_lines, ['Key="Value"'])

    def test_create_key_value_pairs(self):
        # Note: the second test case purposefully has an '=' in the value (the parser only assumes the key has no '=')
        input_lines = ['Key1="Value1"', 'Key2="Value=2"']

        expected_output = [('Key1', 'Value1'), ('Key2', 'Value=2')]
        self.assertEqual(expected_output, self.parser.create_key_value_pairs(input_lines))



    def test_insert_values_into_dict(self):
        original_dict = {'major.minor.A': "init_A_value",
                         'major.minor.B': "init_B_value",
                         'major.minor.C': "init_C_value"
                         }
        key_value_tuples = [('A', '"final_A_value"'), ("B", '"final_B_value"')]

        substituted_dict = self.parser.insert_values_into_dict(original_dict, key_value_tuples)

        # The final dictionary should have new values for A and B, which C's value unchanged
        expected_dict = {'major.minor.A': "final_A_value",
                         'major.minor.B': "final_B_value",
                         'major.minor.C': "init_C_value"
                         }

        self.assertEqual(expected_dict, substituted_dict)


'''
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
'''

if __name__ == '__main__':
    unittest.main()
