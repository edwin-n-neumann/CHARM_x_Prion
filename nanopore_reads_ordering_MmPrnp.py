#!/usr/bin/env python
# coding: utf-8

# ## This script will allow us to re-order nanopore reads based on methylation status of specific CpGs

# ### Edwin Neumann
# ### 12/31/2023

# In[1]:


# import packages
import numpy as np
import pandas as pd
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-e', '--extracted_CpGs')
parser.add_argument('-oo', '--output_order')
parser.add_argument('-i', '--input_bam')
parser.add_argument('-ob', '--output_bam')
args = parser.parse_args()


# In[2]:


# import the modkit extract file with read-specific modified base calls
#extracted_CpGs = pd.read_csv('./experiments/231216_Prnp_KRABless/mouse_980-1/20231216_1354_MN44120_FAY08737_953dd6cb/fast5_pass/pass/extracted.tsv', sep='\t')
extracted_CpGs = pd.read_csv(args.extracted_CpGs, sep='\t')
extracted_CpGs.head()


# In[3]:


# define the window we want to consider (genomic coordinates)
start = 131751500
end = 131752500


# define minimum base quality to consider
qual_min = 9


# In[4]:


# filter extracted_CpGs to just the region of interest and base_qual above threshold
filter_mask = (extracted_CpGs['ref_position'].between(start, end)) & (extracted_CpGs['base_qual'] >= 9)
extracted_CpGs_filtered = extracted_CpGs[filter_mask]
extracted_CpGs_filtered.head()


# In[5]:


# calculate the fraction CpGs methylated for each read
# Group by 'read_id' and perform the calculations
CpGs_calculated = extracted_CpGs_filtered.groupby('read_id').apply(lambda group: pd.Series({
    'fraction_modified': (group['mod_qual'] >= 0.7).mean()
}))

# Reset the index to make 'read_id' a regular column
CpGs_calculated.reset_index(inplace=True)

# now sort the reads by fraction modified, highest to lowest
CpGs_calculated = CpGs_calculated.sort_values(by='fraction_modified', ascending=False).reset_index(drop=True)
CpGs_calculated.head()


# In[6]:


# now save the txt file of the ordered reads
#CpGs_calculated['read_id'].to_csv('./experiments/231216_Prnp_KRABless/mouse_980-1/20231216_1354_MN44120_FAY08737_953dd6cb/fast5_pass/pass/reads_to_keep_ordered.txt', header=False, index=False)
CpGs_calculated['read_id'].to_csv(args.output_order, header=False, index=False)
#print("File 'reads_to_keep_ordered.txt' saved.")


# In[8]:


# rename the reads in the filtered bam file based on the order of reads_to_keep_ordered

def rename_reads(input_bam_filename, output_bam_filename, order_filename):
    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_bam_filename, 'rb') as input_bam:
        # Create a dictionary to map old read names to new names
        read_name_mapping = {}
        
        # Read the order from 'reads_to_keep_ordered.txt'
        with open(order_filename, 'r') as order_file:
            ordered_reads = [line.strip() for line in order_file.readlines()]
        
        # Assign new names starting from 1
        for new_name, old_name in enumerate(ordered_reads, start=100):
            read_name_mapping[old_name] = str(new_name)

        # Open the output BAM file for writing
        with pysam.AlignmentFile(output_bam_filename, 'wb', header=input_bam.header) as output_bam:
            # Iterate through the reads in the input BAM file
            for read in input_bam:
                # Check if the read name is in the mapping
                if read.query_name in read_name_mapping:
                    # Update the read name
                    read.query_name = read_name_mapping[read.query_name]
                    # Write the modified read to the output BAM file
                    output_bam.write(read)

# get file paths 
#order_filename = './experiments/231216_Prnp_KRABless/mouse_980-1/20231216_1354_MN44120_FAY08737_953dd6cb/fast5_pass/pass/reads_to_keep_ordered.txt'
#input_bam_filename = './experiments/231216_Prnp_KRABless/mouse_980-1/20231216_1354_MN44120_FAY08737_953dd6cb/fast5_pass/pass/filtered_output.bam'
#output_bam_filename = './experiments/231216_Prnp_KRABless/mouse_980-1/20231216_1354_MN44120_FAY08737_953dd6cb/fast5_pass/pass/filtered_ordered_output.bam'

# Replace 'input.bam', 'output.bam', and 'reads_to_keep_ordered.txt' with your actual file names
rename_reads(args.input_bam, args.output_bam, args.output_order)


# In[ ]:




