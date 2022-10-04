//
// Created by yonghanyu2 on 8/10/2018.
//

#include "frag_io.h"
#include "type.h"
#include <iostream>
#include "util.h"


template<class ContainerT>
void tokenize(const std::string &str, ContainerT &tokens, const std::string &delimiters, bool trimEmpty)
{
    std::string::size_type pos, lastPos = 0, length = str.length();

    using value_type = typename ContainerT::value_type;
    using size_type  = typename ContainerT::size_type;

    while (lastPos < length + 1)
    {
        pos = str.find_first_of(delimiters, lastPos);
        if (pos == std::string::npos)
        {
            pos = length;
        }

        if (pos != lastPos || !trimEmpty)
            tokens.push_back(value_type(str.data() + lastPos, (size_type) pos - lastPos));

        lastPos = pos + 1;
    }
}

//credit to Marius on stackoverflow:
// https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string


FragmentReader::FragmentReader(const char *file_name)
        : window_start(0), window_end(0), nxt_window_start(0), nxt_window_set(false), intended_window_end(0)
{
    this->prev_chr_var_count = 0;
    this->frag_file.exceptions(std::fstream::failbit | std::fstream::badbit);
    this->buffer.reserve(100);
    try { this->frag_file.open(file_name, std::fstream::in); }
    catch (std::fstream::failure &e)
    {
        std::cerr << "Fail opening fragment file: " << file_name << std::endl
                  << "Check whether the file exits or you have the permission to read." << std::endl;
        exit(-1);
    }
}

FragmentReader::~FragmentReader()
{
    frag_file.close();
}


bool FragmentReader::get_next_pe(Fragment &fragment)
{
    std::fstream::streampos curr_pos;
    try
    {
        curr_pos = this->tell();
        this->frag_file.peek();

        int index_idx = 2;
        int estimated_buffer_len = 2;

        if (NEW_FORMAT)
        {
            index_idx = 5;
            estimated_buffer_len = 5; 
        }

        std::string line;
        this->buffer.clear();

        if (!std::getline(this->frag_file, line))
            return false;

        tokenize(line, this->buffer, " ", true);
        //EOF
        auto token_size = buffer.size();
        if (token_size < 2)
        {
            std::string message = "detected truncated fragment file, exit now.";
            logging(std::cerr, message);
            exit(1);
        }

        int no_blx = std::stoi(this->buffer[0]);

        estimated_buffer_len += no_blx * 2 + 2;

        if (this->buffer.size() < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, the truncated line is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }

        std::string &name = this->buffer[1];
        uint idx_start = std::stol(this->buffer[index_idx]) - 1;

        //new chromosome
        if (idx_start >= curr_chr_var_count + prev_chr_var_count)
        {
            this->seek(curr_pos);
            return false;
        }

        //new phasing window
        if (idx_start >= this->window_end)
        {
            this->seek(this->nxt_window_start);
            return false;
        }

        //set reading position for next phasing window
        if (idx_start >= this->intended_window_end && !this->nxt_window_set)
        {
            this->nxt_window_start = curr_pos;
            this->nxt_window_set = true;
        }

        fragment.read_qual = std::stod(this->buffer.back()) / -10;
        std::string &bs_qual = this->buffer[token_size - 2];

        uint bs_ix = 0;
        uint ix;
        for (int i = 0; i < no_blx; i++)
        {
            ix = std::stol(this->buffer[index_idx + 2 * i]) - 1 - this->prev_chr_var_count; //0-based     //potential problem here
            std::string &blk = this->buffer[2 * i + index_idx + 1];
            for (char &c : blk)
            {
                snp_info t = std::make_pair(ix++, std::make_pair(c - '0', this->cal_base_qual( bs_qual[bs_ix++] )));
                fragment.insert(t);
            }
        }
        fragment.update_start_end();
        fragment.type = FRAG_NORMAL;
        return true;
    }
    catch (const std::ios_base::failure &e)
    {
        return false;
    }
}


bool FragmentReader::get_next_tenx(Fragment &fragment)
{
    std::fstream::streampos curr_pos;
    try
    {
        this->frag_file.peek();
        curr_pos = this->tell();


        std::string line;

        this->buffer.clear();

        if (!std::getline(this->frag_file, line))
            return false;
        tokenize(line, this->buffer, " ", true);
        //EOF
  
        auto token_size = buffer.size();
        if (token_size < 2)
        {
            std::string message = "detected truncated fragment file, exit now.";
            logging(std::cerr, message);
            exit(1);
        }

        int index_idx = 5;
        int estimated_buffer_len = 5; 

        int no_blx = std::stoi(this->buffer[0]);

        estimated_buffer_len += no_blx * 2 + 4;

        if (this->buffer.size() < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, the truncated line is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }

        std::string &name = this->buffer[1];
        uint idx_start = std::stol(this->buffer[index_idx]) - 1;
        

        //new chromosome
        if (idx_start >= curr_chr_var_count + prev_chr_var_count)
        {
            this->seek(curr_pos);
            return false;
        }


        //new phasing window
        if (idx_start >= this->window_end)
        {
            this->seek(this->nxt_window_start);
            return false;
        }

        //set reading position for next phasing window
        if (idx_start >= this->intended_window_end && !this->nxt_window_set)
        {
            this->nxt_window_start = curr_pos;
            this->nxt_window_set = true;
        }

        int rescued = std::stoi(this->buffer[token_size - 2]);
        float dm = std::stof(this->buffer[token_size - 1]);
        fragment.read_qual = std::stod(this->buffer[token_size - 3]) / -10;
        fragment.barcode = this->buffer[3];
        fragment.rescued = rescued == 1;

        fragment.dm = dm;
        std::string &bs_qual = this->buffer[token_size - 4];

        uint bs_ix = 0;
        uint ix;
        for (int i = 0; i < no_blx; i++)
        {
            ix = std::stol(this->buffer[index_idx + 2 * i]) - 1 - this->prev_chr_var_count; //0-based
            std::string &blk = this->buffer[2 * i + index_idx + 1];
            for (char &c : blk)
            {
                snp_info t = std::make_pair(ix++, std::make_pair(c - '0', this->cal_base_qual(bs_qual[bs_ix++])));
                fragment.insert(t);
            }
        }
        fragment.update_start_end();
        fragment.type = FRAG_10X;
        return true;
    }
    catch (const std::ios_base::failure &e)
    {
        return false;
    }
}

// new format, col 2, data type: 1, col 3 mate index, col 4: insertion size, col 5: content
bool FragmentReader::get_next_hic(Fragment &fragment)
{
    std::string line;
    try
    {
        std::fstream::streampos curr_pos = this->tell();
        this->buffer.clear();

        if (!std::getline(this->frag_file, line))
            return false;
        tokenize(line, this->buffer, " ", true);
        //EOF
        auto token_size = buffer.size();
        if (token_size < 2)
        {
            std::string message = "detected truncated fragment file, exit now.";
            logging(std::cerr, message);
            exit(1);
        }

        int index_idx = 5;
        int estimated_buffer_len = 5; 

        int no_blx = std::stoi(this->buffer[0]);

        estimated_buffer_len += no_blx * 2 + 2;

        if (this->buffer.size() < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, the truncated line is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }


        std::string &name = this->buffer[1];
        uint idx_start = std::stol(this->buffer[index_idx]) - 1;
        uint insertion_size = std::stol(this->buffer[4]);
        //new chromosome
        if (idx_start >= curr_chr_var_count + prev_chr_var_count)
        {
            this->seek(curr_pos);
            return false;
        }

        //new phasing window
        if (idx_start >= this->window_end)
        {
            this->seek(this->nxt_window_start);
            return false;
        }

        //set reading position for next phasing window
        if (idx_start >= this->intended_window_end && !this->nxt_window_set)
        {
            this->nxt_window_start = curr_pos;
            this->nxt_window_set = true;
        }

        fragment.read_qual = std::stod(this->buffer.back()) / -10;
        fragment.insertion_size = insertion_size;
        std::string &bs_qual = this->buffer[token_size - 2];

        uint bs_ix = 0;
        uint ix;
        for (int i = 0; i < no_blx; i++)
        {
            ix = std::stol(this->buffer[index_idx + 2 * i]) - 1 - this->prev_chr_var_count; //0-based
            std::string &blk = this->buffer[2 * i + index_idx + 1];
            for (char &c : blk)
            {
                snp_info t = std::make_pair(ix++, std::make_pair(c - '0', this->cal_base_qual(bs_qual[bs_ix++])));
                fragment.insert(t);
            }
        }
        fragment.update_start_end();
        fragment.type = FRAG_HIC;
        return true;
    }
    catch (const std::ios_base::failure &e)
    {
        return false;
    }
}

bool FragmentReader::get_next_nanopore(Fragment &fragment)
{
    std::fstream::streampos curr_pos;
    try
    {
        this->frag_file.peek();
        curr_pos = this->tell();
        int index_idx = 2;
        int estimated_buffer_len = 2;

        if (NEW_FORMAT)
        {
            index_idx = 5;
            estimated_buffer_len = 5; 
        }

        std::string line;

        this->buffer.clear();

        if (!std::getline(this->frag_file, line))
            return false;
        tokenize(line, this->buffer, " ", true);
        //EOF
        auto token_size = buffer.size();
        if (token_size < 2)
        {
            std::string message = "detected truncated fragment file, exit now.";
            logging(std::cerr, message);
            exit(1);
        }

        int no_blx = std::stoi(this->buffer[0]);

        estimated_buffer_len += no_blx * 2 + 2;

        if (this->buffer.size() < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, the truncated line is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }

        std::string &name = this->buffer[1];
        uint idx_start = std::stol(this->buffer[2]) - 1;

        //new chromosome
        if (idx_start >= curr_chr_var_count + prev_chr_var_count)
        {
            this->seek(curr_pos);
            return false;
        }


        //new phasing window
        if (idx_start >= this->window_end)
        {
            this->seek(this->nxt_window_start);
            return false;
        }

        //set reading position for next phasing window
        if (idx_start >= this->intended_window_end && !this->nxt_window_set)
        {
            this->nxt_window_start = curr_pos;
            this->nxt_window_set = true;
        }


        fragment.read_qual = std::stod(this->buffer.back()) / -10;


        std::string &bs_qual = this->buffer[token_size - 2];

        uint bs_ix = 0;
        uint ix;
        for (int i = 0; i < no_blx; i++)
        {
            ix = std::stol(this->buffer[2 + 2 * i]) - 1 - this->prev_chr_var_count; //0-based     //potential problem here
            std::string &blk = this->buffer[2 * i + 3];
            for (char &c : blk)
            {
                snp_info t = std::make_pair(ix++, std::make_pair(c - '0', this->cal_base_qual(bs_qual[bs_ix++])));
                fragment.insert(t);
            }
        }
        fragment.update_start_end();
        fragment.type = FRAG_NORMAL;
        return true;
    }
    catch (const std::ios_base::failure &e)
    {
        return false;
    }
}

bool FragmentReader::get_next_pacbio(Fragment &fragment)
{
    std::fstream::streampos curr_pos;
    try
    {
        this->frag_file.peek();
        curr_pos = this->tell();
        int index_idx = 2;
        int estimated_buffer_len = 2;

        if (NEW_FORMAT)
        {
            index_idx = 5;
            estimated_buffer_len = 5; 
        }

        std::string line;

        this->buffer.clear();

        if (!std::getline(this->frag_file, line))
            return false;
        tokenize(line, this->buffer, " ", true);
        //EOF
        auto token_size = buffer.size();
        
        if (token_size < 2)
        {
            std::string message = "detected truncated fragment file, exit now.";
            logging(std::cerr, message);
            exit(1);
        }

        int no_blx = std::stoi(this->buffer[0]);
        estimated_buffer_len += no_blx * 2 + 2;

        if (this->buffer.size() < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, the truncated line is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }

        std::string &name = this->buffer[1];
        uint idx_start = std::stol(this->buffer[2]) - 1;

        //new chromosome
        if (idx_start >= curr_chr_var_count + prev_chr_var_count)
        {
            this->seek(curr_pos);
            return false;
        }


        //new phasing window
        if (idx_start >= this->window_end)
        {
            this->seek(this->nxt_window_start);
            return false;
        }

        //set reading position for next phasing window
        if (idx_start >= this->intended_window_end && !this->nxt_window_set)
        {
            this->nxt_window_start = curr_pos;
            this->nxt_window_set = true;
        }

        fragment.read_qual = std::stod(this->buffer.back()) / -10;


        std::string &bs_qual = this->buffer[token_size - 2];

        uint bs_ix = 0;
        uint ix;
        for (int i = 0; i < no_blx; i++)
        {
            ix = std::stol(this->buffer[2 + 2 * i]) - 1 - this->prev_chr_var_count; //0-based     //potential problem here
            std::string &blk = this->buffer[2 * i + 3];
            for (char &c : blk)
            {
                snp_info t = std::make_pair(ix++, std::make_pair(c - '0', this->cal_base_qual(bs_qual[bs_ix++])));
                fragment.insert(t);
            }
        }
        fragment.update_start_end();
        fragment.type = FRAG_NORMAL;
        return true;
    }
    catch (const std::ios_base::failure &e)
    {
        return false;
    }
}


bool FragmentReader::get_next_hybrid(Fragment &fragment)
{
    std::fstream::streampos curr_pos;
    try
    {
        curr_pos = this->tell();
        this->frag_file.peek();

        int index_idx = 5;
        int estimated_buffer_len = 5;

        std::string line;
        this->buffer.clear();

        if (!std::getline(this->frag_file, line))
            return false;

        tokenize(line, this->buffer, " ", true);
        //EOF
        auto token_size = buffer.size();
        if (token_size < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, exit now.";
            logging(std::cerr, message);
            exit(1);
        }

        int data_type = std::stoi(this->buffer[2]);

        int no_blx = std::stoi(this->buffer[0]);
        if (data_type == 0)
        {    
            estimated_buffer_len += no_blx * 2 + 2;
            fragment.type = FRAG_NORMAL;
        }
        else if (data_type == 1)  //Hi-C
        {
            estimated_buffer_len += no_blx * 2 + 2;
            fragment.type = FRAG_HIC;
        }
        //TODO: make 10x available in hybrid mode 
        else if (data_type == 2)
        {
            estimated_buffer_len += no_blx * 2 + 4;
            fragment.type = FRAG_10X;
            std::string message = "detected 10x fragment, which is not supported in hybrid mode: " + line ;
            logging(std::cerr, message);
            
            exit(1);
        }
        else 
        {
            std::string message = "detected fragment with unrecognized data type, the line with error is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }

        if (this->buffer.size() < estimated_buffer_len)
        {
            std::string message = "detected truncated fragment file, the truncated line is: " + line ;
            logging(std::cerr, message);
            exit(1);
        }


        std::string &name = this->buffer[1];
        uint idx_start = std::stol(this->buffer[index_idx]) - 1;

        //new chromosome
        if (idx_start >= curr_chr_var_count + prev_chr_var_count)
        {
            this->seek(curr_pos);
            return false;
        }

        //new phasing window
        if (idx_start >= this->window_end)
        {
            this->seek(this->nxt_window_start);
            return false;
        }

        //set reading position for next phasing window
        if (idx_start >= this->intended_window_end && !this->nxt_window_set)
        {
            this->nxt_window_start = curr_pos;
            this->nxt_window_set = true;
        }

        fragment.read_qual = std::stod(this->buffer.back()) / -10;
        std::string &bs_qual = this->buffer[token_size - 2];

        uint bs_ix = 0;
        uint ix;
        for (int i = 0; i < no_blx; i++)
        {
            ix = std::stol(this->buffer[index_idx + 2 * i]) - 1 - this->prev_chr_var_count; //0-based     //potential problem here
            std::string &blk = this->buffer[2 * i + index_idx + 1];
            for (char &c : blk)
            {
                snp_info t = std::make_pair(ix++, std::make_pair(c - '0', this->cal_base_qual( bs_qual[bs_ix++] )));
                fragment.insert(t);
            }
        }
        fragment.update_start_end();
        return true;
    }
    catch (const std::ios_base::failure &e)
    {
        return false;
    }
}


BEDReader::BEDReader(const char *in)
: tbx(nullptr), inf(nullptr), iter(nullptr)
{
    inf = hts_open(in, "r");
    tbx = tbx_index_load(in);
    tmp = {0, 0, nullptr};
}

BEDReader::~BEDReader()
{
    hts_close(inf);
    free(tmp.s);
    tbx_destroy(tbx);
    hts_itr_destroy(iter);
}

int BEDReader::jump_to_region(const char *chromo, uint32_t begin, uint32_t end)
{
    int tid = tbx_name2id(tbx, chromo);
    hts_itr_destroy(iter); iter = nullptr;
    iter = tbx_itr_queryi(tbx, tid, begin, end);
    if (iter == nullptr)
        return -1;
    return 0;
}

int BEDReader::read_region_frag(const char *chromo, uint begin, uint end, RegionFragStats *region_frag_stat)
{
    this->jump_to_region(chromo, begin, end);
    FragStat *buffer = new FragStat();
    while (get_next_record(buffer) == 0)
        region_frag_stat->insert(buffer);

    delete buffer;
    return 0;
}

Range::Range(uint start, uint end)
: start(start), end(end)
{}

bool Range::operator==(const Range &rhs) const
{
    return this->start == rhs.start;
}

std::size_t RangeHasher::operator()(const Range &key) const
{
    return std::hash<uint32_t>()(key.start);
}
