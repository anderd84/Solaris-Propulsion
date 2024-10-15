import sys
import subprocess
import os

cwd = os.getcwd()

print("installing stuff :)")
subprocess.check_call([sys.executable, "-m", "pip", "install", "--editable", cwd])
print("done!")