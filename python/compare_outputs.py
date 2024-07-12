import json

def read_json(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        print(f"Content of {file_path}:")
        print(content)
        return json.loads(content)

def compare_json(python_output, cpp_output):
    py_data = read_json(python_output)
    cpp_data = read_json(cpp_output)

    # Compare sequences
    if py_data['sequences'] != cpp_data['sequences']:
        print("Sequences do not match.")
        return False

    # Compare seed patterns
    py_seed_patterns = [sp['type'] for sp in py_data['seed_patterns']]
    cpp_seed_patterns = [sp['type'] for sp in cpp_data['seed_patterns']]
    if py_seed_patterns != cpp_seed_patterns:
        print("Seed patterns do not match.")
        return False

    # Compare test results
    for py_test, cpp_test in zip(py_data['seed_patterns'], cpp_data['seed_patterns']):
        if py_test['type'] != cpp_test['type']:
            print(f"Seed pattern {py_test['type']} does not match C++ result {cpp_test['type']}.")
            return False
        for py_result, cpp_result in zip(py_test['tests'], cpp_test['tests']):
            if py_result['method'] != cpp_result['method'] or py_result['status'] != cpp_result['status']:
                print(f"Test result for {py_test['type']} does not match. Method: {py_result['method']}.")
                return False

    print("Outputs are identical")
    return True


if __name__ == "__main__":
    python_output_file = 'python_output.n'
    cpp_output_file = 'cpp_output.jn'
    compare_json(python_output_file, cpp_output_file)
