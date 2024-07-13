from memory_profiler import profile
from statistics import mean
from statistics import median
import sys
import glob
from helpers import command_time

@profile
def main():
    if len(sys.argv) != 2:
        exit(f"Usage: {sys.argv[0]} [LANGUAGE]")

    python_base_command = "../python/cplex_env/bin/python ../python/main.py {} {}"
    cpp_base_command = "../cpp/build/ISSH -si {} -q {} -test single"

    if sys.argv[1] == "py":
        base_command = python_base_command
    elif sys.argv[1] == "cpp":
        base_command = cpp_base_command
    else:
        raise Exception("Invalid Language")

    inputs = ["reads_800.fa", "small_test.fa"]
    seeds = glob.glob("../seeds/*")

    for input_file in inputs:
        input_file = "../inputs/" + input_file
        for seed in seeds:
            command_time(base_command.format(input_file, seed))


if __name__ == '__main__':
    main()
