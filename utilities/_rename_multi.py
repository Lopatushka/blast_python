import os
import argparse
'''
    Rename .ab1 files. Test script
'''

def parse_arguments():
    parser = argparse.ArgumentParser(description="Rename files.")
    parser.add_argument('--directory', type=str, help='Directory containing .ab1 files')
    parser.add_argument('--date', type=str, help='Date. Example: "2024-07-31"')
    parser.add_argument('--type', type=str, help='Type of probe. Example: "CS"')
    return parser.parse_args()

def rename_files(directory, date, type):
    for old_name in os.listdir(directory):
        new_name = "Plate-" + date + "_" + type + "_" + old_name
        li = new_name.split("_")
        probe_name = "-".join([li[2], li[3]])
        li[2:4] = [probe_name]
        new_name = "_".join(li)
        os.rename(
                os.path.join(directory, old_name),
                os.path.join(directory, new_name)
            )
        print(f"File renamed from '{old_name}' to '{new_name}' in directory '{directory}'.")

def main():
    args = parse_arguments()
    dir = args.directory
    date = args.date
    type = args.type

    rename_files(directory = dir, date = date, type = type)

if __name__ == "__main__":
    main()


