import sys
import unittest

from Bashwrapper import *

class BashWrappertest(unittest.TestCase):

    def test_exception(self):
        with self.assertRaises(Exception):
            Bash("awk")

    def test_stdout(self):
        cmd = Bash("echo \'hello world\'")
        self.assertEqual(cmd.stdout(), 'hello world\n')

    def test_stderr(self):
        with self.assertRaises(Exception):
            cmd = Bash("awk")
            self.assertTrue(cmd.stderr())

if __name__ == '__main__':
    unittest.main()
