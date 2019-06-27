import pyBigWig

MAX_CHROM_SIZE = 250000000

files = [
    'random_test.bw',
    'sequential_test_ascending.bw',
    'sequential_test_same.bw',
    'skip_test_same.bw'
]

for file in files:
    print(f"Testing {file}")

    bw = pyBigWig.open(file)

    for i in range(100, 10000, 500):
        print(bw.stats('chr1', 1, i))

    print()
