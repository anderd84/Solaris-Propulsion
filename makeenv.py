import sys
import subprocess
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

print("installing stuff :)")
subprocess.check_call([sys.executable, "-m", "pip", "install", "--editable", dir_path])
print("done!")