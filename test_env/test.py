import sys
import re

def main():
    path_to_summary = sys.argv[1]
    scaffold="NZ_CP048008.1"    
    
    with open (path_to_summary, 'r') as file:
        file_content = file.read()
        regex = f"({scaffold}\n[-]*\n)(Cluster \\d*)"
        print(regex)
        result = re.search(regex, file_content)
        print(result.group(2))


if __name__ == "__main__":
    main()
