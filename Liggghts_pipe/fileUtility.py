import sys
import os
import glob

def my_copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def get_filename_without_extension(filepath):
  """
  From Gemini 2.0
  Extracts the filename from a given filepath, removing the extension.

  Args:
    filepath: The full path to the file.

  Returns:
    The filename without the extension, or the original filename if no extension is found.
  """
  filename = os.path.basename(filepath)  # Get the filename with extension
  filename_without_extension, _ = os.path.splitext(filename) #Split on the extension
  return filename_without_extension