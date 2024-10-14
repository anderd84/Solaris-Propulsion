import sys
import subprocess

print("installing stuff :)")
subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
print("done!")