//
// Created by yyh on 3/4/2019.
//

#include "phaser.h"
#include "htslib/vcf.h"
#include "type.h"
#include <iostream>
#include <cstdlib>
#include <unordered_map>

// TODO clarify between variant count and block count

Phaser::Phaser(std::vector<option::Option> &options)
{
    threshold = pow(10, -7);
    if (options[TENX])
        op = op_mode::TENX;
    else if (options[HIC])
        op = op_mode::HIC;
    else if (options[PACBIO])
        op = op_mode::PACBIO;
    else if (options[NANOPORE])
        op = op_mode::NANOPORE;
    else
        op = op_mode::PE;
    //op_mode frag_file_op_mode = detect_frag_file_type(options[FRAGMENT].arg);


    frvcf = new VCFReader(options[VCF].arg);
    fwvcf = new VCFWriter(frvcf->header, options[OUT].arg);
    frfrag = new FragmentReader(options[FRAGMENT].arg);

    if (options[WINDOW_SIZE].arg == nullptr)
        window_size = 200;
    else
        window_size = uint(atoi(options[WINDOW_SIZE].arg));
    if (options[WINDOW_OVERLAP].arg == nullptr)
        overlap = 60;
    else
        overlap = int(atoi(options[WINDOW_OVERLAP].arg));

    coverage = 30;  //deprecated
    max_barcode_spanning_length = 60000;    
    recursive_limit = 15;       //deprecated

    if (options[TENX])
        frbed = new BEDReader(options[STATS].arg);

    bool use_secondary = false;

    spectral = new Spectral(frfrag, frbed, op, threshold, coverage, max_barcode_spanning_length, use_secondary);
}

Phaser::~Phaser()
{
    delete frvcf;
    delete fwvcf;
    delete frfrag;
    delete spectral;
    delete frbed;
}


int Phaser::load_contig_records(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    while (true)
    {
        ptr_ResultforSingleVariant result = std::make_shared<ResultforSingleVariant>();
        if (this->frvcf->get_next_record_contig(result) != 0)
            break;
        chromo_phaser->results_for_variant.push_back(result);
    }


    chromo_phaser->variant_count = chromo_phaser->results_for_variant.size();
    for (int i = 0; i < chromo_phaser->variant_count; i++)
    {
        chromo_phaser->variant_to_block_id[i] = i;
    }
    chromo_phaser->block_count = chromo_phaser->variant_count;
    return status;
}

int Phaser::load_contig_blocks(ChromoPhaser *chromo_phaser)
{
    int status = 0;
    this->load_contig_records(chromo_phaser);
    
    std::unordered_map<uint, uint> ps2block_ids;
    uint block_count = 0;
    for (int i = 0; i < chromo_phaser->variant_count; i++)
    {
        auto result = chromo_phaser->results_for_variant[i];
        uint ps = result->ps;
        if (ps == 0) //not phased 
        {
            chromo_phaser->variant_to_block_id[i] = block_count;
            block_count++;
        }
        else {      //phased
            if (ps2block_ids.count(ps) == 0)
            {   // not met before
                ps2block_ids[ps] = block_count;
                chromo_phaser->variant_to_block_id[i] = ps2block_ids[ps];
                block_count++;
            }
            else {
                chromo_phaser->variant_to_block_id[i] = ps2block_ids[ps];
            }
        }
    }
    chromo_phaser->block_count = block_count;
    return status;
}

void Phaser::phasing()
{
    uint prev_variant_count = 0;

    for (uint rid = 0; rid < frvcf->contigs_count; rid++)
    {
        if (frvcf->jump_to_contig(rid) != 0)
            break;
        ChromoPhaser *chromo_phaser = new ChromoPhaser(rid, frvcf->contigs[rid], overlap, window_size);
        load_contig_records(chromo_phaser);
        chromo_phaser->construct_phasing_window_initialize();
        frfrag->set_prev_chr_var_count(prev_variant_count);
        spectral->set_chromo_phaser(chromo_phaser);
        phasing_by_chrom(chromo_phaser->variant_count, chromo_phaser);

        fwvcf->write_nxt_contigs(frvcf->contigs[rid].data(), chromo_phaser, *frvcf);
        //write vcf
        prev_variant_count += chromo_phaser->variant_count;
        spectral->release_chromo_phaser();
        delete chromo_phaser;
    }
}

void Phaser::phasing_by_chrom(uint var_count, ChromoPhaser *chromo_phaser)
{
    frfrag->set_curr_chr_var_count(var_count);

    while (chromo_phaser->phased->rest_blk_count > 0)
    {
        if (chromo_phaser->phased->rest_blk_count > chromo_phaser->variant_count)
            break;
        if (op == op_mode::TENX)
            chromo_phaser->phased->update_phasing_info(max_barcode_spanning_length);
        else
            chromo_phaser->phased->update_phasing_info();
        spectral->solver();
        //std::cout << chromo_phaser->phased->phased_blk_count << std::endl;

    }

    if (op == op_mode::HIC)
        spectral->hic_poss_solver();

}

void Phaser::phase_HiC_recursive(ChromoPhaser *chromo_phaser)
{
    uint prev_blk_count;
    for(int j = 0; j < recursive_limit; j++)
    {
        prev_blk_count = chromo_phaser->get_block_count();
        chromo_phaser->construct_phasing_window_r_initialize();
        while (chromo_phaser->phased->rest_blk_count > 0)
        {
            chromo_phaser->phased->update_phasing_info();
            spectral->clean();
            spectral->solver_recursive();
        }
        if (prev_blk_count == chromo_phaser->get_block_count())
            break;
        prev_blk_count = chromo_phaser->get_block_count();
    }
}

