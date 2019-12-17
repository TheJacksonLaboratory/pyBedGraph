from pyBedGraph import BedGraph
import os
import csv

FOLDER_LOC = '/media/hirow/extra/jax/data/pybedgraph'

stats = []

for folder in os.listdir(FOLDER_LOC):
    for filename in os.listdir(f'{FOLDER_LOC}/{folder}'):
        print(folder, filename)
        bedgraph = BedGraph(f'/media/hirow/extra/jax/data/chrom_sizes/{folder}.chrom.sizes', f'{FOLDER_LOC}/{folder}/{filename}', 'chr1')
        sample_name = filename.split('.')[0]
        sample = {}
        sample['name'] = sample_name
        chrom = bedgraph.chromosome_map['chr1']

        sample['total_coverage'] = chrom.total_coverage
        sample['num_samples'] = chrom.num_samples
        sample['avg_chrom_value'] = chrom.avg_chrom_value
        sample['avg_interval_value'] = chrom.avg_interval_value
        sample['avg_interval_size'] = chrom.avg_interval_size
        sample['num_intervals'] = chrom.num_intervals

        stats.append(sample)

csv_columns = list(stats[0].keys())

with open('bedgraph_stats.csv', 'w') as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=csv_columns)
    writer.writeheader()
    for data in stats:
        writer.writerow(data)
