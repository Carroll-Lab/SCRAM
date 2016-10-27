import unittest
import scram_modules.analysis_helper as ah

class TestAHMethods(unittest.TestCase):

    def test_single_file_output(self):
        in_file = "/home/foo/bar/my_file.fa"
        out_name= "my_file"
        self.assertEqual(ah.single_file_output(in_file), out_name)

if __name__ == '__main__':
    unittest.main()