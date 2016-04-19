if __name__ == "__main__":
    import sys
    import subprocess
    import shlex
    from operator import ior
    sys.exit(reduce(ior, [subprocess.check_call(shlex.split("bash pipeline_with_qc.sh %s %s" % (arg, arg.replace("R1", "R2"))))
        for arg in sys.argv[1:]]))
