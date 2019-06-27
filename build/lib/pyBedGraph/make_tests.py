import random

NUMB_TEST_CASES = 10000
MAX_CHROM_SIZE = 250000000

with open('data/test_files/sequential_test_same.bedgraph', 'w') as test:
    output = ""
    value = 0.1
    for i in range(NUMB_TEST_CASES):
        start = i * 100 + 1
        end = start + 99
        output += f"chr1\t{start}\t{end}\t{value}\n"

    test.write(output)

with open('data/test_files/sequential_test_ascending.bedgraph', 'w') as test:
    output = ""
    value = 0.1
    for i in range(NUMB_TEST_CASES):
        start = i * 100 + 1
        end = start + 99
        output += f"chr1\t{start}\t{end}\t{round(value * i, 1)}\n"

    test.write(output)

with open('data/test_files/skip_test_same.bedgraph', 'w') as test:
    output = ""
    value = 0.1
    for i in range(1, NUMB_TEST_CASES * 2, 200):
        start = i * 100 + 1
        end = start + 99
        output += f"chr1\t{start}\t{end}\t{value}\n"

    test.write(output)

with open('data/test_files/random_test.bedgraph', 'w') as test:
    output = ""
    current_index = 1
    for i in range(NUMB_TEST_CASES):
        start = current_index
        end = start + random.randint(50, 2000)
        current_index = end + random.randint(50, 2000)
        value = random.uniform(0.1, 100.0)
        output += f"chr1\t{start}\t{end}\t{value}\n"

        if end > MAX_CHROM_SIZE:
            print(f"Cannot create a chromosome above size: {MAX_CHROM_SIZE}")
            exit(-1)

    test.write(output)

with open('data/test_files/myChrom.sizes', 'w') as test:
    test.write(f"chr1\t{MAX_CHROM_SIZE}")

with open('intervals.txt', 'w') as test:
    output = ""
    for i in range(NUMB_TEST_CASES):
        start = i * 100 + 1
        end = start + 99
        output += f"{start}\t{end}\n"

    test.write(output)

