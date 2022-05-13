"""Test for the core.py program."""

import unittest

class TestCore(unittest.TestCase):
    """Test methods class."""

    def test(self):
        """Square test."""
        self.assertAlmostEqual(2.**2, 4.)

if __name__ == "__main__":
    unittest.main()
