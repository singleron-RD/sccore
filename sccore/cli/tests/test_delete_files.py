import unittest
import os
import time
from sccore.cli.delete_files import should_delete_file, delete_files


class TestDeleteFiles(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory structure with test files
        self.test_dir = "test_directory"
        os.makedirs(self.test_dir, exist_ok=True)

        # Create test files with different modification times
        self.old_file = os.path.join(self.test_dir, "old_file.bam")
        with open(self.old_file, "w") as f:
            f.write("This is an old file.")
        # Set modification time to 100 days ago
        old_time = time.time() - (100 * 24 * 60 * 60)
        os.utime(self.old_file, (old_time, old_time))

        self.recent_file = os.path.join(self.test_dir, "recent_file.bam")
        with open(self.recent_file, "w") as f:
            f.write("This is a recent file.")
        # Set modification time to 10 days ago
        recent_time = time.time() - (10 * 24 * 60 * 60)
        os.utime(self.recent_file, (recent_time, recent_time))

        # File with keyword in name
        self.keyword_file = os.path.join(self.test_dir, "test.bam")
        with open(self.keyword_file, "w") as f:
            f.write("This file should match the keyword criteria.")
        # Set modification time to 100 days ago
        old_time = time.time() - (100 * 24 * 60 * 60)
        os.utime(self.keyword_file, (old_time, old_time))

    def tearDown(self):
        # Remove the test directory and files after each test
        for root, dirs, files in os.walk(self.test_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(self.test_dir)

    def test_should_delete_old_file(self):
        # Test that old files are marked for deletion
        self.assertTrue(should_delete_file(self.old_file))

    def test_should_not_delete_recent_file(self):
        # Test that recent files are not marked for deletion
        self.assertFalse(should_delete_file(self.recent_file))

    def test_should_delete_file_with_keyword(self):
        # Test that files with keywords are marked for deletion
        self.assertTrue(should_delete_file(self.keyword_file))

    def test_delete_files(self):
        # Test the delete_files function
        delete_files([self.test_dir])
        # Old and keyword files should be deleted, recent file should still exist
        self.assertFalse(os.path.exists(self.old_file))
        self.assertFalse(os.path.exists(self.keyword_file))
        self.assertTrue(os.path.exists(self.recent_file))


if __name__ == "__main__":
    unittest.main()
