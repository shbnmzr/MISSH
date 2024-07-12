import json
import sys


def is_executable():
    if len(sys.argv) != 3:
        exit(f"Usage: {sys.argv[0]} [FASTA] [SEED]")


def log_output(filename: str, data: dict) -> None:
    with open(filename, 'w') as output_file:
        json.dump(data, output_file, indent=4)
