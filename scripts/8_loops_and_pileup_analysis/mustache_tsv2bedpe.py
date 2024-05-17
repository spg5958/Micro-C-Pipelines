import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Convert mustache .tsv to .bedpe')
parser.add_argument('-infile', required=True, help='mustache tsv file')
parser.add_argument('-outfile', required=True, help='output file')
args = parser.parse_args()

infile_path = args.infile
outfile_path = args.outfile

print(infile_path)
print(outfile_path)

df=pd.read_csv(infile_path, sep="\t")
print(df)

out_df=df.iloc[:,:6]
print(out_df)
out_df.to_csv(outfile_path, sep="\t", header=None, index=False)
