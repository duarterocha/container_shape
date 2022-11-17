from pyoomph.expressions.units import *
from pyoomph.utils.paramscan import *

if __name__ == "__main__":
    # Create a parameter scanner, give the script to run and the max number of processes to run simultaneously
    scanner = ParallelParameterScan("container_shape.py", max_procs=6)

    for gamma in numpy.logspace(-1,1,20):
        sim = scanner.new_sim("Gamma_" + str(gamma))
        sim.aspect_ratio = gamma

    # Run all (and rerun also already finished sims)
    scanner.run_all(skip_done=True)

    shell_file = '/Users/duarterocha17/Desktop/container_shape/merge_output_files.sh'
    P = subprocess.Popen(shell_file)