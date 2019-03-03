# Usage: use python2
# Command: python2 script.py
#
# install XlsxWriter by pip install XlsxWriter or anaconda

from os import listdir, getcwd
from os.path import join, dirname, realpath
import sys
import pandas as pd

script_path = dirname(realpath(__file__))
project_name = '16p_11.2'
tsv_output = 'output.tsv'

headers = ['Project',
           'Sample_name',
           'Cutadapt_Total_read_pairs_processed',
           'Cutadapt_Read1_with_adapter_count',
           'Cutadapt_Read1_with_adapter_percent',
           'Cutadapt_Read2_with_adapter_count',
           'Cutadapt_Read2_with_adapter_percent',
           'Cutadapt_Pairs_too_short_count',
           'Cutadapt_Pairs_too_short_percent',
           'Cutadapt_Pairs_passing_filters_count',
           'Cutadapt_Pairs_passing_filters_percent',
           'Cutadapt_Total_basepairs_processed',
           'Cutadapt_Read1_basepairs_processed',
           'Cutadapt_Read2_basepairs_processed',
           'Cutadapt_Total_quality_trimmed_basepairs_count',
           'Cutadapt_Total_quality_trimmed_basepairs_percent',
           'Cutadapt_Read1_quality_trimmed_basepairs',
           'Cutadapt_Read2_quality_trimmed_basepairs',
           'Cutadapt_Total_written_basepairs_count',
           'Cutadapt_Total_written_basepairs_percent',
           'Cutadapt_Read1_written_basepairs',
           'Cutadapt_Read2_written_basepairs',
           'Cutadapt_First_read_adapter',
           'Cutadapt_First_read_type',
           'Cutadapt_First_read_length',
           'Cutadapt_First_read_trimmed_times',
           'Cutadapt_First_read_bases_preceding_removed_adapters_percent_A',
           'Cutadapt_First_read_bases_preceding_removed_adapters_percent_C',
           'Cutadapt_First_read_bases_preceding_removed_adapters_percent_G',
           'Cutadapt_First_read_bases_preceding_removed_adapters_percent_T',
           'Cutadapt_First_read_bases_preceding_removed_adapters_percent_Other',
           'Cutadapt_Second_read_adapter',
           'Cutadapt_Second_read_type',
           'Cutadapt_Second_read_length',
           'Cutadapt_Second_read_trimmed_times',
           'Cutadapt_Second_read_bases_preceding_removed_adapters_percent_A',
           'Cutadapt_Second_read_bases_preceding_removed_adapters_percent_C',
           'Cutadapt_Second_read_bases_preceding_removed_adapters_percent_G',
           'Cutadapt_Second_read_bases_preceding_removed_adapters_percent_T',
           'Cutadapt_Second_read_bases_preceding_removed_adapters_percent_' +
           'Other',
           'STAR_uniquely_mapped_percent',
           'STAR_num_splices',
           'STAR_num_GCAG_splices',
           'STAR_insertion_length',
           'STAR_deletion_length',
           'STAR_unmapped_tooshort_percent',
           'STAR_avg_mapped_read_length',
           'STAR_deletion_rate',
           'STAR_mismatch_rate',
           'STAR_avg_input_read_length',
           'STAR_num_ATAC_splices',
           'STAR_num_annotated_splices',
           'STAR_num_GTAG_splices',
           'STAR_uniquely_mapped',
           'STAR_multimapped_toomany',
           'STAR_unmapped_mismatches',
           'STAR_unmapped_mismatches_percent',
           'STAR_total_reads',
           'STAR_unmapped_other',
           'STAR_insertion_rate',
           'STAR_unmapped_other_percent',
           'STAR_multimapped_percent',
           'STAR_multimapped',
           'STAR_num_noncanonical_splices',
           'STAR_unmapped_tooshort',
           'STAR_multimapped_toomany_percent',
           'Picard_MarkDuplicates_Unmapped_reads',
           'Picard_MarkDuplicates_Estimated_library_size',
           'Picard_MarkDuplicates_Read_pair_optical_duplicates',
           'Picard_MarkDuplicates_Read_pairs_examined',
           'Picard_MarkDuplicates_Unpaired_reads_examined',
           'Picard_MarkDuplicates_Percent_duplication',
           'Picard_MarkDuplicates_Read_pair_duplicates',
           'Picard_MarkDuplicates_Unpaired_read_duplicates',
           'Picard_CollectRnaSeqMetrics_Coding_bases',
           'Picard_CollectRnaSeqMetrics_Percent_coding_bases',
           'Picard_CollectRnaSeqMetrics_Intronic_bases',
           'Picard_CollectRnaSeqMetrics_Percent_intronic_bases',
           'Picard_CollectRnaSeqMetrics_Intergenic_bases',
           'Picard_CollectRnaSeqMetrics_Percent_intergenic_bases',
           'Picard_CollectRnaSeqMetrics_Ribosomal_bases',
           'Picard_CollectRnaSeqMetrics_Percent_ribosomal_bases',
           'Picard_CollectRnaSeqMetrics_Percent_mRNA_bases',
           'Picard_CollectRnaSeqMetrics_UTR_bases',
           'Picard_CollectRnaSeqMetrics_Percent_UTR_bases',
           'Picard_CollectRnaSeqMetrics_Correct_strand_reads',
           'Picard_CollectRnaSeqMetrics_Percent_correct_strand_reads',
           'Picard_CollectRnaSeqMetrics_Incorrect_strand_reads',
           'Picard_CollectRnaSeqMetrics_Ignored_reads',
           'Picard_CollectRnaSeqMetrics_Percent_usable_bases',
           'Picard_CollectRnaSeqMetrics_Median_3prime_bias',
           'Picard_CollectRnaSeqMetrics_Median_5prime_bias',
           'Picard_CollectRnaSeqMetrics_Median_5prime_to_3prime_bias',
           'Picard_CollectRnaSeqMetrics_Median_cv_coverage',
           'Picard_CollectRnaSeqMetrics_Pass_filter_bases',
           'Picard_CollectRnaSeqMetrics_Pass_filter_aligned_bases',
           'Picard_CollectInsertSizeMetrics_Median_insert_size',
           'Picard_CollectInsertSizeMetrics_Median_absolute_deviation',
           'Picard_CollectInsertSizeMetrics_Min_insert_size',
           'Picard_CollectInsertSizeMetrics_Max_insert_size',
           'Picard_CollectInsertSizeMetrics_Mean_insert_size',
           'Picard_CollectInsertSizeMetrics_Standard_deviation',
           'Picard_CollectInsertSizeMetrics_Read_pairs',
           'Picard_CollectInsertSizeMetrics_Pair_orientation',
           'Picard_CollectInsertSizeMetrics_Width_of_10_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_20_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_30_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_40_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_50_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_60_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_70_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_80_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_90_percent',
           'Picard_CollectInsertSizeMetrics_Width_of_99_percent',
           'Picard_CollectAlignmentSummaryMetrics_Total_reads',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_reads',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_reads_aligned',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_aligned_bases',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_high_quality_' +
           'aligned_reads',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_high_quality_' +
           'aligned_bases',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_high_quality_' +
           'aligned_q20_bases',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_high_quality_' +
           'error_rate',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_high_quality_' +
           'median_mismatches',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_mismatch_rate',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_indel_rate',
           'Picard_CollectAlignmentSummaryMetrics_Pass_filter_noise_reads',
           'Picard_CollectAlignmentSummaryMetrics_Percent_reads_aligned_in_' +
           'pairs',
           'Picard_CollectAlignmentSummaryMetrics_Percent_pass_filter_reads',
           'Picard_CollectAlignmentSummaryMetrics_Percent_pass_filter_reads_' +
           'aligned',
           'Picard_CollectAlignmentSummaryMetrics_Percent_chimeras',
           'Picard_CollectAlignmentSummaryMetrics_Percent_adapter',
           'Picard_CollectAlignmentSummaryMetrics_Reads_aligned_in_pairs',
           'Picard_CollectAlignmentSummaryMetrics_Mean_read_length',
           'Picard_CollectAlignmentSummaryMetrics_Bad_cycles',
           'Picard_CollectAlignmentSummaryMetrics_Strand_balance',
           'Picard_CollectAlignmentSummaryMetrics_Category',
           'Picard_CollectGcBiasMetrics_Total_clusters',
           'Picard_CollectGcBiasMetrics_Aligned_reads',
           'Picard_CollectGcBiasMetrics_AT_dropout',
           'Picard_CollectGcBiasMetrics_GC_dropout',
           'RNA-SeQC_End1_mapping_rate',
           'RNA-SeQC_End1_mismatch_rate',
           'RNA-SeQC_End1_sense',
           'RNA-SeQC_End1_percent_sense',
           'RNA-SeQC_End1_antisense',
           'RNA-SeQC_End2_mapping_rate',
           'RNA-SeQC_End2_mismatch_rate',
           'RNA-SeQC_End2_sense',
           'RNA-SeQC_End2_percent_sense',
           'RNA-SeQC_End2_antisense',
           'RNA-SeQC_Mapping_rate',
           'RNA-SeQC_Base_mismatch_rate',
           'RNA-SeQC_Intragenic_rate',
           'RNA-SeQC_Intergenic_rate',
           'RNA-SeQC_Exonic_rate',
           'RNA-SeQC_Intronic_rate',
           'RNA-SeQC_rRNA_rate',
           'RNA-SeQC_Mapped_unique_rate_of_total',
           'RNA-SeQC_Unique_rate_of_mapped',
           'RNA-SeQC_Mapped_unique',
           'RNA-SeQC_Duplication_rate_of_mapped',
           'RNA-SeQC_3prime_norm',
           'RNA-SeQC_5prime_norm',
           'RNA-SeQC_Num_covered_5prime',
           'RNA-SeQC_Cumulative_gap_length',
           'RNA-SeQC_Num_gaps',
           'RNA-SeQC_Gap_percent',
           'RNA-SeQC_Genes_detected',
           'RNA-SeQC_Transcripts_detected',
           'RNA-SeQC_Mapped',
           'RNA-SeQC_Mapped_pairs',
           'RNA-SeQC_Unpaired_reads',
           'RNA-SeQC_Total_purity_filtered_reads_sequenced',
           'RNA-SeQC_Estimated_library_size',
           'RNA-SeQC_Failed_vendor_QC_check',
           'RNA-SeQC_Chimeric_pairs',
           'RNA-SeQC_Split_reads',
           'RNA-SeQC_Read_length',
           'RNA-SeQC_Expression_profiling_efficiency',
           'RNA-SeQC_Alternative_alignments',
           'RNA-SeQC_Fragment_length_standard_deviation',
           'RNA-SeQC_Fragment_length_mean',
           'RNA-SeQC_Mean_per_base_coverage',
           'RNA-SeQC_Mean_cv',
           'RNA-SeQC_rRNA',
           'featureCounts_Unassigned_ambiguity',
           'featureCounts_Unassigned_mapping_quality',
           'featureCounts_Percent_assigned',
           'featureCounts_Unassigned_nonjunction',
           'featureCounts_Unassigned_duplicate',
           'featureCounts_Unassigned_chimera',
           'featureCounts_Unassigned_unmapped',
           'featureCounts_Assigned',
           'featureCounts_Unassigned_multimapping',
           'featureCounts_Unassigned_secondary',
           'featureCounts_Unassigned_no_features',
           'featureCounts_Unassigned_fragment_length',
           'featureCounts_Total',
           'ChrM_percent']


def combine():
    sample_names = []
    sample_dirs = []
    if 'samples.txt' not in listdir(script_path):
        print('Error: missing samples.txt')
        exit()
    with open(join(script_path, 'samples.txt'), 'r') as f:
        lines = f.read().splitlines()
    lines = sorted(lines)
    for line in lines:
        sample_names.append(line.split()[0])
        sample_dirs.append(line.split()[1])
    table = [[title] for title in headers]
    for index in range(len(sample_names)):
        table[0].append(project_name)
        table[1].append(sample_names[index])
    table = add_cutadapt(table, sample_names, sample_dirs)
    table = add_star(table, sample_names, sample_dirs)
    table = add_picard(table, sample_names, sample_dirs)
    table = add_rnaseqc(table, sample_names, sample_dirs)
    table = add_fc(table, sample_names, sample_dirs)
    table = add_chrm(table, sample_names, sample_dirs)
    table = pd.DataFrame(list(zip(*table)))
    table.to_csv(join(script_path, tsv_output), sep='\t', index=False,
                 header=False)


def extract_percent(raw):
    for i in range(len(raw)):
        if raw[i] == '%':
            number = raw[1:i]
            if number == '0.0':
                return '0'
            return number


def add_cutadapt(table, names, dirs):
    files = []
    for i in range(len(names)):
        file_name = names[i] + '.Cutadapt.out'
        if file_name not in listdir(dirs[i]):
            print('Error: missing cutadapt files: ' + names[i])
            exit()
        files.append(join(dirs[i], file_name))
    for file_index in range(len(names)):
        with open(files[file_index], 'r') as f:
            while True:
                if "=== Summary ===" in f.readline():
                    break
            f.readline()
            lines = f.read().splitlines()
            index = 1
            for line in lines:
                if line == '':
                    continue
                if line[:3] == '===':
                    break
                if line.endswith(')'):
                    index += 1
                    if line.split()[-2] == 'bp':
                        number = line.split()[-3]
                        table[index].append(''.join(number.split(',')))
                    else:
                        number = line.split()[-2]
                        table[index].append(''.join(number.split(',')))
                    index += 1
                    table[index].append(extract_percent(line.split()[-1]))
                elif line.endswith('bp'):
                    index += 1
                    number = line.split()[-2]
                    table[index].append(''.join(number.split(',')))
                else:
                    index += 1
                    number = line.split()[-1]
                    table[index].append(''.join(number.split(',')))
        for i in range(22, 40):
            table[i].append('NA')
    return table


def add_star(table, names, dirs):
    files = []
    for i in range(len(names)):
        file_name = names[i] + '.STARLog.final.out'
        if file_name not in listdir(dirs[i]):
            print('Error: missing STAR files: ' + name)
            exit()
        files.append(join(dirs[i], file_name))
    for file_index in range(len(names)):
        with open(files[file_index], 'r') as f:
            lines = f.read().splitlines()
        for line_index in range(len(lines)):
            if lines[line_index] == '' or lines[line_index].endswith(':'):
                continue
            lines[line_index] = lines[line_index].split('\t')[1].strip()
        table[40].append(lines[9][:-1])
        table[41].append(lines[11])
        table[42].append(lines[14])
        table[43].append(lines[21])
        table[44].append(lines[19])
        too_short_percent = lines[29][:-1]
        table[45].append(too_short_percent)
        table[46].append(lines[10])
        table[47].append(lines[18][:-1])
        table[48].append(lines[17][:-1])
        table[49].append(lines[6])
        table[50].append(lines[15])
        table[51].append(lines[12])
        table[52].append(lines[13])
        table[53].append(lines[8])
        table[54].append(lines[25])
        mismatch_percent = lines[28][:-1]
        table[56].append(mismatch_percent)
        total_reads = lines[5]
        table[57].append(total_reads)
        table[59].append(lines[20][:-1])
        other_percent = lines[30][:-1]
        table[60].append(other_percent)
        table[61].append(lines[24][:-1])
        table[62].append(lines[23])
        table[63].append(lines[16])
        table[65].append(lines[26][:-1])
        table[55].append(round(float(mismatch_percent) * float(total_reads)
                               * 0.01))
        table[58].append(round(float(other_percent) * float(total_reads) *
                               0.01))
        table[64].append(round(float(too_short_percent) *
                               float(total_reads) * 0.01))
    return table


def add_picard(table, names, dirs):
    md_files = []
    is_files = []
    rna_files = []
    as_files = []
    gc_files = []
    for i in range(len(names)):
        folder_name = names[i] + '.Picard'
        files = listdir(join(dirs[i], folder_name))
        md_file = names[i] + '.marked_dup_metrics.txt'
        is_file = names[i] + '.insert_size_metrics.txt'
        rna_file = names[i] + '.RNA_Metrics.txt'
        as_file = names[i] + '.CollectalignmentSummaryMetricsoutput.txt'
        gc_file = names[i] + '.summary_metrics.txt'
        count = 0
        for file in files:
            if file == md_file:
                md_files.append(join(dirs[i], folder_name, md_file))
            	count += 1
            elif file == is_file:
                is_files.append(join(dirs[i], folder_name, is_file))
            	count += 1
            elif file == rna_file:
                rna_files.append(join(dirs[i], folder_name, rna_file))
            	count += 1
            elif file == as_file:
                as_files.append(join(dirs[i], folder_name, as_file))
            	count += 1
            elif file == gc_file:
                gc_files.append(join(dirs[i], folder_name, gc_file))
            	count += 1
        if count != 5:
            print('Error: missing Picard files: ' + name)
            exit()
    for file_index in range(len(names)):
        with open(md_files[file_index], 'r') as f:
            lines = f.read().splitlines()
        titles = lines[6].split('\t')
        numbers = lines[7].split('\t')
        data = dict(zip(titles, numbers))
        for index in range(66, 74):
            temp = headers[index].split('_')[2:]
            if temp[0] == 'Percent':
                key = '_'.join(temp).upper()
                number = data[key]
                percent = float(number) * 100
                table[index].append(str(percent))
                continue
            key = '_'.join(temp).upper()
            table[index].append(data[key])
        with open(rna_files[file_index], 'r') as f:
            lines = f.read().splitlines()
        titles = lines[6].split('\t')
        numbers = lines[7].split('\t')
        data = dict(zip(titles, numbers))
        for index in range(74, 96):
            temp = headers[index].split('_')[2:]
            if temp[0] == 'Percent':
                temp[0] = 'PCT'
                key = '_'.join(temp).upper()
                number = data[key]
                percent = float(number) * 100
                table[index].append(str(percent))
                continue
            if temp[0] == 'Pass':
                temp[0] = 'PF'
                del temp[1]
            key = '_'.join(temp).upper()
            table[index].append(data[key])
        with open(is_files[file_index], 'r') as f:
            lines = f.read().splitlines()
        titles = lines[6].split('\t')
        numbers = lines[7].split('\t')
        data = dict(zip(titles, numbers))
        for index in range(96, 114):
            temp = headers[index].split('_')[2:]
            key = '_'.join(temp).upper()
            table[index].append(data[key])
        with open(as_files[file_index], 'r') as f:
            lines = f.read().splitlines()
        titles = lines[6].split('\t')
        numbers = lines[9].split('\t')
        data = dict(zip(titles, numbers))
        for index in range(114, 136):
            temp = headers[index].split('_')[2:]
            if temp[0] == 'Percent':
                temp[0] = 'PCT'
                if temp[1] == 'pass':
                    temp[1] = 'PF'
                    del temp[2]
                key = '_'.join(temp).upper()
                number = data[key]
                percent = float(number) * 100
                table[index].append(str(percent))
                continue
            if temp[0] == 'Pass':
                temp[0] = 'PF'
                del temp[1]
                if temp[1] == 'high':
                    temp[1] = 'HQ'
                    del temp[2]
            key = '_'.join(temp).upper()
            table[index].append(data[key])
        with open(gc_files[file_index], 'r') as f:
            lines = f.read().splitlines()
        titles = lines[6].split('\t')
        numbers = lines[7].split('\t')
        data = dict(zip(titles, numbers))
        for index in range(136, 140):
            temp = headers[index].split('_')[2:]
            key = '_'.join(temp).upper()
            table[index].append(data[key])
    return table


def add_rnaseqc(table, names, dirs):
    files = []
    for i in range(len(names)):
        folder1 = names[i] + '.RNA-SeQC'
        folder2 = names[i] + '.testReport'
        if 'metrics.tsv' not in listdir(join(dirs[i], folder1, folder2)):
            print('Error: missing RNA_SeQC files: ' + names[i])
            exit()
        files.append(join(dirs[i], folder1, folder2, 'metrics.tsv'))
    for file_index in range(len(names)):
        with open(files[file_index], 'r') as f:
            lines = f.read().splitlines()
        titles = lines[0].split('\t')
        numbers = lines[1].split('\t')
        data = dict(zip(titles, numbers))
        for index in range(140, 185):
            temp = headers[index].split('_')[1:]
            if temp[0] == '3prime' or temp[0] == '5prime':
                temp[0] = temp[0].replace('prime', '\'')
            if 'End' in temp[0]:
                temp[0] = 'End ' + temp[0][-1]
            if temp[0] == 'Num':
                if temp[1] == 'covered':
                    temp[0] = 'No.'
                    temp[2] = '5\''
                else:
                    temp[0] = 'Num.'
            if temp[0] == 'Failed':
                temp[2] = 'QC'
            temp[0] = temp[0].replace('Cumulative', 'Cumul.')
            if temp[-1] == 'coverage':
                temp[-1] = 'Cov.'
            if temp[-1] == 'deviation':
                temp[2] = 'StdDev'
                del temp[3]
            if temp[-1] == 'cv':
                temp[-1] = 'CV'
            if temp[0] == 'Alternative':
                temp[1] = 'Aligments'
            for i in range(len(temp)):
                temp[i] = temp[i].replace('percent', '%')
                if temp[i] == 'rRNA':
                    break
                elif temp[i] != 'of':
                    temp[i] = temp[i][0].upper() + temp[i][1:]
            key = ' '.join(temp)
            table[index].append(data[key])
    return table


def add_fc(table, names, dirs):
    files = []
    for i in range(len(names)):
        folder_name = names[i] + '.FeatureCounts'
        file_name = names[i] + '.counts.txt.summary'
        if file_name not in listdir(join(dirs[i], folder_name)):
            print('Error: missing featureCounts files: ' + names[i])
            exit()
        files.append(join(dirs[i], folder_name, file_name))
    for file_index in range(len(names)):
        with open(files[file_index], 'r') as f:
            lines = f.read().splitlines()
        ambiquity = lines[12].split('\t')[1]
        table[185].append(ambiquity)
        mappingquiality = lines[3].split('\t')[1]
        table[186].append(mappingquiality)
        nonjunction = lines[9].split('\t')[1]
        table[188].append(nonjunction)
        duplicate = lines[6].split('\t')[1]
        table[189].append(duplicate)
        chimera = lines[4].split('\t')[1]
        table[190].append(chimera)
        unmapped = lines[2].split('\t')[1]
        table[191].append(unmapped)
        assigned = lines[1].split('\t')[1]
        table[192].append(assigned)
        multimapping = lines[7].split('\t')[1]
        table[193].append(multimapping)
        secondary = lines[8].split('\t')[1]
        table[194].append(secondary)
        nofeatures = lines[10].split('\t')[1]
        table[195].append(nofeatures)
        fragmentlength = lines[5].split('\t')[1]
        table[196].append(fragmentlength)
        total = int(ambiquity) + int(mappingquiality) + int(nonjunction) + \
            int(duplicate) + int(chimera) + int(unmapped) + int(assigned) + \
            int(multimapping) + int(secondary) + int(nofeatures) + \
            int(fragmentlength)
        table[197].append(str(total))
        table[187].append(str(int(assigned) / total * 100))
    return table


def add_chrm(table, names, dirs):
    files = []
    for i in range(len(names)):
        folder_name = names[i] + '.RNA-SeQC'
        file_name = names[i] + '.Mt.txt'
        if file_name not in listdir(join(dirs[i], folder_name)):
            print('Error: missing Mt files: ' + names[i])
            exit()
        files.append(join(dirs[i], folder_name, file_name))
    for file_index in range(len(names)):
        with open(files[file_index], 'r') as f:
            lines = f.read().splitlines()
        number = lines[0]
        if number[0] == '.':
            number = '0' + number
        table[198].append(number)
    return table


if __name__ == '__main__':
    combine()
