import os
import argparse
import glob
from colorama import init, Fore, Style
import logging
from typing import List, Tuple, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

init(autoreset=True)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def format_size(size_in_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(size_in_bytes)
    unit_index = 0

    while size >= 1024 and unit_index < len(units) - 1:
        size /= 1024
        unit_index += 1

    return f"{size:.2f} {units[unit_index]}"


def get_file_size(file_path: str) -> Optional[Tuple[str, int]]:
    try:
        if os.path.isfile(file_path):
            size = os.path.getsize(file_path)
            return (file_path, size)
        else:
            logging.warning(f"{file_path}: Not a file or doesn't exist")
            return None
    except Exception as e:
        logging.error(f"Error processing {file_path}: {str(e)}")
        return None


def get_file_sizes(path: str, pattern: Optional[str] = None) -> List[Tuple[str, int]]:
    file_info = []

    try:
        if os.path.isfile(path):
            with open(path, "r") as f:
                files = [line.strip() for line in f]
        elif os.path.isdir(path) or pattern:
            if pattern:
                files = glob.glob(os.path.join(path, pattern))
            else:
                files = [
                    os.path.join(root, f)
                    for root, _, filenames in os.walk(path)
                    for f in filenames
                ]
        else:
            logging.error(f"{path}: Not a file or directory")
            return []

        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(get_file_size, file_path) for file_path in files]
            for future in as_completed(futures):
                result = future.result()
                if result:
                    file_info.append(result)

    except Exception as e:
        logging.error(f"Error processing path {path}: {str(e)}")

    return sorted(file_info, key=lambda x: x[1], reverse=True)


def print_file_sizes(file_info: List[Tuple[str, int]]) -> None:
    for path, size in file_info:
        formatted_size = format_size(size)

        if size < 1e9:  # Less than 1GB
            color = Fore.RED
        elif size < 10e9:  # Less than 10GB
            color = Fore.YELLOW
        else:
            color = Fore.RESET

        print(f"{path}: {color}{formatted_size}{Style.RESET_ALL}")


def main():
    parser = argparse.ArgumentParser(
        description="Get file sizes from a list of paths in a file, from a directory, or matching a pattern."
    )
    parser.add_argument(
        "path",
        help="Path to the file containing the list of file paths, a directory, or a pattern",
    )
    parser.add_argument(
        "-p", "--pattern", help="File pattern to match (e.g., '*.bam')", default=None
    )
    args = parser.parse_args()

    file_info = get_file_sizes(args.path, args.pattern)
    print_file_sizes(file_info)


if __name__ == "__main__":
    main()
