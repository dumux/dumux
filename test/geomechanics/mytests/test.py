#!/usr/bin/env python3
import os
import subprocess
testList = ["test_mohr_circle",
            "test_point_line_distance",
            "test_stress_drop_law_params",
            "test_stress_drop_law",
            "test_constantcementmodel"]

dir = os.path.dirname(__file__)
os.chdir(dir)

for test in testList:
    if os.path.exists(test):
        os.remove(test)
    subprocess.run(["make", f"{test}"], check=True)


for test in testList:
    subprocess.run([f"./{test}"])
    print(f"\n test {test} passed. \n")