import psutil
import threading
import os
import time
import subprocess


def command_in_thread(command: str):
    def thread_function():
        os.system(command)

    thread = threading.Thread(target=thread_function)
    thread.start()
    thread.join()


def command_time(command: str):
    start_time = time.time()
    # os.system(command)
    run_command(command)
    return round((time.time() - start_time) * 1000, 2)


def run_command(command: str):
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def memory():
    process = psutil.Process()
    return process.memory_info().rss
