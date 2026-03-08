import gzip
import pandas as pd

def create_dbsnp_mapping(vcf_path, output_tsv):
    """
    Extract rsid, chrom, pos, ref, alt from dbSNP VCF.
    vcf_path: path to .vcf.gz file
    output_tsv: path to output tab-separated file
    """
    data = []
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            # dbSNP VCF may have multiple alts; we'll take the first
            alt = alt.split(',')[0]
            data.append([rsid, chrom, pos, ref, alt])
    
    df = pd.DataFrame(data, columns=['rsid', 'chrom', 'pos', 'ref', 'alt'])
    df['pos'] = df['pos'].astype(int)
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Created mapping file with {len(df)} entries")

# Usage
create_dbsnp_mapping('GCF_000001405.40.gz', 'dbsnp_b38_mapping.tsv')