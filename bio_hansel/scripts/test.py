import io
import os
import pandas as pd



def read_vcf():
    with open('core.vcf', 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_table(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str}
    ).rename(columns={'#CHROM': 'CHROM'})

def main():
    data_frame=read_vcf()
    # print(data_frame)
    data_frame=data_frame[data_frame['ALT'].str.len()>1]
    print(data_frame)

    header = """##fileformat=VCFv4.1
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##source=myImputationProgramV3.1
    ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele.">
    #CHROM POS ID REF ALT QUAL FILTER INFO
    """

    output_VCF = "myfile.vcf"
    with open(output_VCF, 'w') as vcf:
      vcf.write(header)
    
    data_frame.to_csv(output_VCF, sep="\t", mode='a', index=False)


if __name__ == '__main__':
    main()