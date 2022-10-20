import time
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--use-caching", action="store_true")
args = vars(parser.parse_args())

print("Will run the grid geometry test once with the new and the old implementation")
print("Monitor your system's memory consumption to compare the memory footprint of the two")

exe_prefix = "with" if args["use_caching"] else "no"
exe_new = f"test_grid_geometry_performance_{exe_prefix}_cache_new_implementation"
exe_old = f"test_grid_geometry_performance_{exe_prefix}_cache_old_implementation"

subprocess.run(["make", exe_new], check=True)
subprocess.run(["make", exe_old], check=True)
print("Sleeping 5s for your cpus to cool down after compiling/before execution")
time.sleep(5)

subprocess.run([f"./{exe_new}"], check=True)
print("Sleeping 5s for you to check your memory usage")
time.sleep(5)
subprocess.run([f"./{exe_old}"], check=True)
