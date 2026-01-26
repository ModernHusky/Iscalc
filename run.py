import sys

import integral.state

def process_file(filename):
    try:
        with open(filename, 'r') as file:
            content = file.read()
        # Call the check_actions function with the file content
        result = integral.state.check_actions(content, print_lines=True)
        return result
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run.py <filename>")
    else:
        filename = sys.argv[1]
        process_file(filename)
