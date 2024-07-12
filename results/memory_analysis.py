from helpers import memory
from helpers import command_in_thread

initial_memory = memory()
command_in_thread("../cplex_env/bin/python ../main.py")
during_memory = memory()
