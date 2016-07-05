from collections import defaultdict as dd


def get_single_doc(path):
    for line in open(path):
       if "Mean Singleton Coverage" in line:
           return float(line.split()[-1])
    raise ValueError("Missing field!")


def get_raw_doc(path):
    for line in open(path):
       if "Mean Raw Coverage" in line:
           return float(line.split()[-1])
    raise ValueError("Missing field!")


def get_mean_doc(path):
    sys.stderr.write("Path: %s.\n" % path)
    for line in open(path):
        print(line)
        if "Mean Collapsed Coverage" in line:
            return float(line.split()[-1])
    raise ValueError("Missing field!")


def get_mean_egfr(path):
    c = s = 0
    for line in open(path):
        if line[0] == "#": continue
        if line.split()[3].startswith("EGFR"):
            c += 1
            s += float(line.split()[4].split(":")[1])
    try: ret = s / c
    except ZeroDivisionError:
        print(path); raise
    return ret

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python %s <bedpath> <bampath> <bampath2> ....\n"
                         "Emits mean doc and egfr to stdout.\n" % sys.argv[0])
        sys.exit(1)
    sys.stdout.write("Name\tMean DOC\tMean Raw DOC\tMean Singleton DOC\tMean EGFR DOC\n")
    for path in sys.argv[1:]:
        try: m, r, s, e = get_mean_doc(path), get_raw_doc(path), get_single_doc(path), get_mean_egfr(path)
        except ValueError:
            print("Failed at ", path)
            raise
        sys.stdout.write("%s\t%f\t%f\t%f\t%f\n" % (path, m, r, s, e))
