//
// Created by yyh on 3/4/2019.
//

#ifndef SPECHAP_PHASER_H
#define SPECHAP_PHASER_H

#include "type.h"
#include "vcf_io.h"
#include "frag_io.h"
#include "optionparser.h"

enum optionIndex
{
    UNKNOWN, HELP, VCF, FRAGMENT, OUT, TENX, HIC, WINDOW_SIZE, COVERAGE, RECURSIVE_LIMIT, NANOPORE, PACBIO, NOSORT,
    MAX_BARCODE_SPANNING_LENGTH, WINDOW_OVERLAP, STATS, NEWFORMAT, USESECONDARY, KEEP_PHASING_INFO
};

class Phaser
{
public:
    Phaser() = default;
    explicit Phaser(std::vector<option::Option> &options);
    ~Phaser();
    void phasing();

private:
    op_mode detect_frag_file_type(std::string file_name);
    void sort_frag_file(std::string file_name, op_mode opMode);
    double threshold;
    VCFReader *frvcf;
    VCFWriter *fwvcf;
    FragmentReader *frfrag;
    BEDReader *frbed;
    uint window_size;
    int overlap;
    Spectral *spectral;
    op_mode op;
    uint recursive_limit;
    int coverage;
    int max_barcode_spanning_length;
    void phasing_by_chrom(uint var_count, ChromoPhaser *chromo_phaser);
    void phase_HiC_recursive(ChromoPhaser *chromo_phaser, std::set<uint> &connected_comp);
    void phase_HiC_poss(ChromoPhaser *chromo_phaser);
    void update_HiC_phasing_window();
    int load_contig_records(ChromoPhaser *chromo_phaser);
    int load_contig_blocks(ChromoPhaser *chromo_phaser);
    bool use_input_ps;
};

#endif //SPECHAP_PHASER_H
