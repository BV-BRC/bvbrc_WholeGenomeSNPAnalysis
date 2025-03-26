import sys
import re

def parse_optimum_k(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r'The optimum value of k is (\d+)', line)
            if match:
                print(match.group(1))
                return
    print("Optimum value of k not found")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
    else:
        parse_optimum_k(sys.argv[1])
