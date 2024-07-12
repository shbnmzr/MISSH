import sys
import glob
from statistics import mean
from statistics import median
from memory_profiler import profile

@profile
def check():
    a = [i for i in range(100000)]
    b = [i * 2 for i in a]
    b = 100 * b
    return b

if __name__ == '__main__':
    check()

def main():
    if len(sys.argv) != 2:
        exit(f"Usage: {sys.argv[0]} [LANGUAGE]")

    python_base_command = "../python/cplex_env/bin/python ../python/main.py {} {}"
    cpp_base_command = "../cpp/build/ISSH -si {} -q {}"

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
            experiments = []
            for i in range(50):
                experiments.append()
            print(f"{sys.argv[1]} {input_file} {seed}: {round(mean(experiments), 2)} {round(median(experiments), 2)}")
