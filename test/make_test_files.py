import random

NUMB_TEST_CASES = 5
MAX_CHROM_SIZE = 30  # mouse chr1

'''with open('test_files/sequential_test_same.bedgraph', 'w') as test:
    output = ""
    value = 0.1
    for i in range(NUMB_TEST_CASES):
        start = i * 100 + 1
        end = start + 99
        output += f"chr1\t{start}\t{end}\t{value}\n"

    test.write(output)

with open('test_files/sequential_test_ascending.bedgraph', 'w') as test:
    output = ""
    value = 0.1
    for i in range(NUMB_TEST_CASES):
        start = i * 100 + 1
        end = start + 99
        output += f"chr1\t{start}\t{end}\t{round(value * i, 1)}\n"

    test.write(output)

with open('test_files/skip_test_same.bedgraph', 'w') as test:
    output = ""
    value = 0.1
    for i in range(1, NUMB_TEST_CASES * 2, 200):
        start = i * 100 + 1
        end = start + 99
        output += f"chr1\t{start}\t{end}\t{value}\n"

    test.write(output)'''

with open('test_files/random_test.bedGraph', 'w') as test:
    output = ""
    current_index = 0
    for i in range(NUMB_TEST_CASES):
        start = current_index
        end = start + random.randint(1, 5)
        current_index = end + random.randint(0, 5)
        value = random.random()
        value = round(value, 1)
        output += f"chr1\t{start}\t{end}\t{value}\n"

        if end > MAX_CHROM_SIZE:
            print(f"Cannot create a chromosome above size: {MAX_CHROM_SIZE}")
            exit(-1)

    test.write(output)

with open('test_files/myChrom.sizes', 'w+') as test:
    test.write(f"chr1\t{MAX_CHROM_SIZE}")

with open('test_files/test_intervals.txt', 'w+') as test:
    output = ""
    for i in range(NUMB_TEST_CASES):
        interval_size = random.randint(1, 5)
        start = random.randint(0, MAX_CHROM_SIZE - interval_size)
        end = start + interval_size
        output += f"chr1\t{start}\t{end}\n"

    test.write(output)
