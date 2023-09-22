"""
This script is written to provide a class for standardizing the output of YACHT. Things like the CAMI profiling format, BIOM format, something compatible with GraphPlAn, etc.
"""

# import libraries
import os, sys
import pandas as pd
import pytaxonkit
import numpy as np
import biom
import argparse
from biom.util import biom_open
from utils import get_cami_profile
from collections import OrderedDict
from loguru import logger
logger.remove()
logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} - {level} - {message}", level="INFO")

class StandardizeYachtOutput:
    """
    Standardize the output of YACHT to a format (options: CAMI, BIOM, GraphPlAn)
    """
    def __init__(self):
        ## set allowable_rank
        self.allowable_rank = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
        
    def __to_cami(self, sample_name):
        """
        Convert the YACHT output to the CAMI format.
        params:sample_name: The sample name shown in header of the CAMI file
        return:None
        """
        ## Find the ncbi lineage
        result = pytaxonkit.lineage(list(set(self.genome_to_taxid['genome_id'].tolist())))
        metadata_df = self.genome_to_taxid.merge(result[['TaxID','Rank','FullLineageTaxIDs','FullLineage','FullLineageRanks']], left_on='taxid', right_on='TaxID').drop(columns=['taxid']).reset_index(drop=True)
        metadata_df['FullLineageTaxIDs'] = metadata_df['FullLineageTaxIDs'].str.replace(';','|')
        metadata_df['FullLineage'] = metadata_df['FullLineage'].str.replace(';','|')
        metadata_df['FullLineageRanks'] = metadata_df['FullLineageRanks'].str.replace(';','|')

        ## select the organisms that YACHT considers to present in the sample
        yacht_res_df = self.yacht_output.copy()
        organism_name_list = yacht_res_df['organism_name'].tolist()
        organism_id_list = [x.split(' ')[0] for x in organism_name_list]

        if len(organism_id_list) == 0:
            logger.error("No organism is detected by YACHT.")
            exit(1)

        ## Merge the YACHT results with the metadata
        selected_organism_metadata_df = metadata_df.query("genome_id in @organism_id_list").reset_index(drop=True)

        ## Summarize the results
        summary_dict = {}
        for row in selected_organism_metadata_df.to_numpy():
            select_index = [index for index, x in enumerate(row[-1].split('|')) if x in self.allowable_rank]
            taxid_list = list(np.array(row[3].split('|'))[select_index])
            lineage_list = list(np.array(row[4].split('|'))[select_index])
            rank_list = list(np.array(row[5].split('|'))[select_index])
            current_lineage = ''
            current_taxid = ''
            for index, (taxid, rank, lineage) in enumerate(zip(taxid_list, rank_list, lineage_list)):
                if index == 0:
                    current_lineage = lineage
                    current_taxid = taxid
                else: 
                    current_lineage = current_lineage + '|' + lineage
                    current_taxid = current_taxid + '|' + taxid
                if taxid not in summary_dict:
                    summary_dict[taxid] = {'RANK':rank, 'TAXPATH':current_taxid, 'TAXPATHSN':current_lineage, 'count':1}
                else:
                    summary_dict[taxid]['count'] += 1

        # calculate percentage
        for taxid in summary_dict:
            summary_dict[taxid]['PERCENTAGE'] = summary_dict[taxid]['count']/len(selected_organism_metadata_df)*100

        ## sort by rank in allowable rank list
        summary_df = pd.DataFrame(summary_dict).T.reset_index().rename(columns={'index':'TAXID'})
        res_df = [summary_df.query(f'RANK == "{rank}"') for rank in self.allowable_rank]
        res_df = pd.concat(res_df).drop(columns=['count']).reset_index(drop=True)
        res_df.columns = ['@@TAXID', 'RANK', 'TAXPATH', 'TAXPATHSN', 'PERCENTAGE']

        ## output summary results
        if len(res_df) == 0:
            logger.error("Unable to convert to a CAMI format. No organism is detected by YACHT.")
            exit(1)
        else:
            output_line_list = []
            output_line_list.append(f"@SampleID:{sample_name}")
            output_line_list.append('@Version:0.9.1')
            output_line_list.append(f"@Ranks:{'|'.join(list(res_df['RANK'].unique()))}")
            output_line_list.append('')
            output_line_list.append('@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE')
            for row in res_df.to_numpy():
                output_line_list.append(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}")
                
        return output_line_list

    def __to_biom(self, sample_name):
        """
        Convert the YACHT output to the BIOM format.
        params:sample_name: The sample name shown in header of the CAMI file
        return:None
        """
        # call the _to_cami function
        try:
            cami_content = self.__to_cami(sample_name)
        except Exception:
            logger.error("Unable to convert to a BIOM format. No organism is detected by YACHT.")
            exit(1)
        samples_list = get_cami_profile(cami_content)
        
        ## convert to BIOM format
        tax_ids_row_mapping = OrderedDict()
        sample_ids = []
        sample_metadata = []
        data = {}
        observ_metadata = {}
        
        for i, sample in enumerate(samples_list):
            sample_id, header, profile = sample

            if sample_id == '':
                sample_id = f'S{i}'
                header['SAMPLEID'] = sample_id

            sample_ids.append(sample_id)
            sample_metadata.append(header)

            for prediction in profile:
                if prediction.taxid not in tax_ids_row_mapping:
                    tax_ids_row_mapping[prediction.taxid] = len(tax_ids_row_mapping)
                data[(tax_ids_row_mapping[prediction.taxid], i)] = prediction.percentage
                observ_metadata[prediction.taxid] = prediction.get_metadata()

        # sort metadata
        observ_metadata = [observ_metadata[key] for key in tax_ids_row_mapping.keys()]

        table = biom.Table(data=data,
                        observation_ids=list(tax_ids_row_mapping.keys()),
                        sample_ids=sample_ids,
                        observation_metadata=observ_metadata,
                        sample_metadata=sample_metadata,
                        table_id='table_0',
                        type='OTU table')

        return table
    
    def __to_graphplan(self, sample_name):
        """
        Convert the YACHT output to the input of GraphPlAn.
        params:sample_name: The sample name shown in header of the CAMI file
        return:None
        """
        
        def tree_to_newick(node):
            if not node:
                return ""
            children = ",".join([f"{child}{tree_to_newick(node[child])}" for child in node])
            return f"({children})"
        
        # call the _to_cami function
        try:
            cami_content = self.__to_cami(sample_name)
        except Exception:
            logger.error("Unable to convert to a GraphPlAn format. No organism is detected by YACHT.")
            exit(1)

        ## convert to GraphPlAn format
        taxa_tree = {}
        annotations = []
        for line in cami_content:
            # skip the header
            if line.startswith("@@") or line.startswith("@") or line == '':
                continue
            
            # skip the invalid line
            if line.startswith("#"):
                ranks = line.split(':')[1].split('|')
                break

            parts = line.split('\t')
            taxpathsn = parts[3].split('|')
            percentage = parts[4]

            # Build the tree structure
            parent = taxa_tree
            for taxon in taxpathsn:
                if taxon not in parent:
                    parent[taxon] = {}
                parent = parent[taxon]

        return tree_to_newick(taxa_tree)
    
    def __savefile(self, output_format, path_to_outdir, content, fileprefix='result'):
        """
        Save the output file.
        params:output_format: The format of the output file. Options: cami, biom, graphplan.
        params:path_to_outdir: The path to the output directory.
        params:content: The content of the output file.
        params:fileprefix: The prefix of the output file.
        """
    
        if output_format == 'cami':
            with open(os.path.join(path_to_outdir, f"{fileprefix}.cami"), 'w') as f:
                for line in content:
                    f.write(f"{line}\n")
        elif output_format == 'biom':
            with biom_open(os.path.join(path_to_outdir, f"{fileprefix}.biom"), 'w') as f:
                content.to_hdf5(f, 'YACHT')
        elif output_format == 'graphplan':
            with open(os.path.join(path_to_outdir, f"{fileprefix}.graphplan.newick"), 'w') as f:
                f.write(content + ';')
        else:
            logger.error(f"{output_format} is not a valid output format. Please choose from cami, biom, graphplan.")
            return

    
    def run(self, yacht_output, genome_to_taxid, path_to_outdir, fileprefix='result', output_format='cami', sample_name='Sample1'):
        """
        The main function to run the format standardization.        
        params:yacht_output: A dataframe containing the YACHT output that only contains the organisms that are considered to be present in the sample.
        params:genome_to_taxid: A dataframe with two columns: genome ID (genome_id) and its corresponding taxid (taxid).
        params:path_to_outdir: The path to the output directory.
        params:fileprefix: The prefix of the output file.
        params:output_format: The format of the output file. Options: cami, biom, graphplan.
        params:sample_name: The sample name shown in header of the output file.
        """
        self.yacht_output = yacht_output
        self.genome_to_taxid = genome_to_taxid
        
        if output_format == 'cami':
            result = self.__to_cami(sample_name)
            logger.info(f"Saving the CAMI format to {os.path.join(path_to_outdir, f'{fileprefix}.cami')}")
        elif output_format == 'biom':
            result = self.__to_biom(sample_name)
            logger.info(f"Saving the BIOM format to {os.path.join(path_to_outdir, f'{fileprefix}.biom')}")
        elif output_format == 'graphplan':
            result = self.__to_graphplan(sample_name)
            logger.info(f"Saving the GraphPlAn tree file (newick format) to {path_to_outdir}.")
        else:
            logger.error(f"{output_format} is not a valid output format. Please choose from cami, biom, graphplan.")
            return
        
        # save the file
        self.__savefile(output_format, path_to_outdir, result, fileprefix)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script convert YACHT output to a format (options: CAMI, BIOM, GraphPlAn).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--yacht_output', type=str, help='Path to the YACHT output excel file.', required=True)
    parser.add_argument('--sheet_name', type=str, help='The sheet name of the YACHT output excel file.', required=True)
    parser.add_argument('--genome_to_taxid', type=str, help='Path to the genome to taxid file. This file is a TSV file \
                        with two columns: genome ID (genome_id) and its corresponding taxid (taxid).', required=True)
    parser.add_argument('--mode', type=str, help='The output format. Options: cami, biom, graphplan, all. Default: cami.', required=False, default='cami')
    parser.add_argument('--sample_name', type=str, help='The sample name shown in header of the file. Default: Sample1.', required=False, default='Sample1')
    parser.add_argument('--outfile_prefix', type=str, help='The prefix of the output file. Default: result.', required=False, default='result')
    parser.add_argument('--outdir', type=str, help='The path to the output directory.', required=True)
    
    # parse the arguments
    args = parser.parse_args()
    yacht_output = args.yacht_output
    sheet_name = args.sheet_name
    genome_to_taxid = args.genome_to_taxid
    mode = args.mode
    sample_name = args.sample_name
    outfile_prefix = args.outfile_prefix
    outdir = args.outdir
    
    # check if the yacht output file exists
    if not os.path.exists(yacht_output):
        logger.error(f"{yacht_output} does not exist.")
        raise ValueError
    
    # check if the genome to taxid file exists
    if not os.path.exists(genome_to_taxid):
        logger.error(f"{genome_to_taxid} does not exist.")
        raise ValueError
    
    # check if the output directory exists and create it if not
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    # load the yacht output
    yacht_output_df = pd.read_excel(yacht_output, sheet_name=sheet_name, engine='openpyxl')
    # converet the first column to string
    yacht_output_df['organism_name'] = yacht_output_df['organism_name'].astype(str)
    
    # load the genome to taxid file
    genome_to_taxid_df = pd.read_csv(genome_to_taxid, sep='\t', header=0)
    # converet the first column to string
    genome_to_taxid_df['genome_id'] = genome_to_taxid_df['genome_id'].astype(str)
    
    # run the standardization
    standardize_yacht_output = StandardizeYachtOutput()
    if mode == 'all':
        standardize_yacht_output.run(yacht_output_df, genome_to_taxid_df, outdir, outfile_prefix, 'cami', sample_name)
        standardize_yacht_output.run(yacht_output_df, genome_to_taxid_df, outdir, outfile_prefix, 'biom', sample_name)
        standardize_yacht_output.run(yacht_output_df, genome_to_taxid_df, outdir, outfile_prefix, 'graphplan', sample_name)
    elif mode == 'cami':
        standardize_yacht_output.run(yacht_output_df, genome_to_taxid_df, outdir, outfile_prefix, 'cami', sample_name)
    elif mode == 'biom':
        standardize_yacht_output.run(yacht_output_df, genome_to_taxid_df, outdir, outfile_prefix, 'biom', sample_name)
    elif mode == 'graphplan':
        standardize_yacht_output.run(yacht_output_df, genome_to_taxid_df, outdir, outfile_prefix, 'graphplan', sample_name)
    else:
        logger.error(f"{mode} is not a valid output format. Please choose from cami, biom, graphplan, all.")
        exit(1)