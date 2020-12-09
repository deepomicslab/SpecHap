#include "bamread.h"

char INT_CIGAROP[] = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};

extern int DATA_TYPE;

int QVoffset = 33;

void print_read_debug(struct alignedread* read)
{
	fprintf(stderr,"%s %d %d IS: %d %d\t",read->readid,read->tid,read->position,read->IS,read->cigs);
	int i=0; for (i=0;i<read->cigs;i++) fprintf(stderr,"%d%c",read->cigarlist[i] >> 4, INT_CIGAROP[read->cigarlist[i]&0xf]); 
	fprintf(stderr,"\n");
}

int fetch_func(const bam1_t *b, void *data, struct alignedread* read) {
    samfile_t *fp = (samfile_t*) data;
    uint32_t *cigar = bam1_cigar(b);
    const bam1_core_t *c = &b->core;
    int i, op, ol;
    read->cigs = 0;
    read->alignedbases = 0;
    read->clipped = 0;
    read->span = 0;
    read->gapped = 0;
    read->cflag = 0;
    read->readlength = b->core.l_qseq;
    read->sequence = (char*) malloc(b->core.l_qseq + 1);
    read->quality = (char*) malloc(b->core.l_qseq + 1);

    uint8_t* sequence = bam1_seq(b);
    uint8_t* quality = bam1_qual(b);
    for (i = 0; i < b->core.l_qseq; i++) read->sequence[i] = bam_nt16_rev_table[bam1_seqi(sequence, i)];
    read->sequence[i] = '\0';
    if (quality[0] == 255) // quality string is missing, 01/29/2014, quality is set to minimum quality value specified using --minq
    {
        for (i = 0; i < b->core.l_qseq; i++) read->quality[i] = (char) (MINQ + 33);
        read->quality[i] = '\0';
    } else {
        for (i = 0; i < b->core.l_qseq; i++) read->quality[i] = (char) (quality[i] + 33);
        read->quality[i] = '\0';
    }
    //fprintf(stderr,"quality |%d| \n",quality[1]);

    read->flag = c->flag;
    read->mquality = c->qual;
    read->position = c->pos + 1;
    read->mateposition = c->mpos + 1;
    read->IS = c->isize;
    read->strand = '+';
    if ((read->flag & 16) == 16) read->strand = '-'; // fixed sept 29 2011

    read->cigarlist = (int*) malloc(sizeof (int)*c->n_cigar);
    //read->cigs = c->n_cigar; 
    int j=0,pc=0;
    for (i = 0; i < c->n_cigar; ++i) {
        op = cigar[i]&0xf;
        ol = cigar[i] >> 4;
        if (op == BAM_CEQUAL || op == BAM_CDIFF) // convert to 'M' cigar format 
	{
		if (pc ==0) read->cigarlist[j++] = ol<<4;   // new 'M' operation
		else read->cigarlist[j-1] += ol<<4; // add to previous cigar 'M' operation 
		pc =1; 
	} 
	else { read->cigarlist[j++] = cigar[i]; pc =0; } 

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            read->alignedbases += ol;
            read->span += ol;
        } else if (op == BAM_CDEL) {
            read->gapped += 1;
            read->span += ol;
        } else if (op == BAM_CINS) {
            read->alignedbases += ol;
            read->gapped += 1;
        } else if (op == BAM_CREF_SKIP) read->span += ol;
        else if (op == BAM_CSOFT_CLIP) read->clipped += ol;
        else if (op == BAM_CHARD_CLIP) {
        } else read->cflag = 1;
    }
    read->cigs = j;
    char* barcode = NULL;
    read->barcode = NULL;
	if ( !(read->flag & 4) && DATA_TYPE == 2) // 10X reads
	{
		barcode = (char *) bam_aux_get(b,"BX");
		if (barcode && barcode != NULL){
            read->barcode = (char*)malloc(strlen(barcode)); //else read->barcode = NULL;
        	for (op=1;op<strlen(barcode);op++){
                read->barcode[op-1] = barcode[op];
            }
            read->barcode[op-1]= '\0';
        }else read->barcode = NULL;
		char *rescue = (char *) bam_aux_get(b, "OM");
		char *diverg = (char *) bam_aux_get(b, "DM");
		read->rescued = false;
		if (rescue != NULL)
        {
		    int32_t r = bam_aux2i(rescue);
            if (read->mquality >= 30)
            {
                if (r < 30)
                    read->rescued = true;
            }
        } else read->rescued = false;
		if (diverg != NULL)
        {
		    char *temp = (char*) malloc(strlen(diverg));
		    int k = 0;
		    for (k = 1; k < strlen(diverg); k++)
            {
		        temp[k-1] = diverg[k];
            }
            temp[k - 1] = '\0';
		    read->dm = atof(temp);
		    free(temp);
        }else read->dm = 0;
	}

    //if (read->mquality >= 60) read->mquality = 60; // cap it at 60 april 18 2012
    read->readid = (char*) malloc(c->l_qname + 1);
    unsigned char* qs = b->data;
    for (i = 0; i < c->l_qname; i++) read->readid[i] = qs[i];
    read->readid[i] = '\0';

    if (c->tid >= 0) read->chrom = fp->header->target_name[c->tid];
    else read->chrom = NULL;
    if (c->mtid >= 0) read->matechrom = fp->header->target_name[c->mtid];
    else read->matechrom = NULL;
    read->tid = c->tid;
    read->mtid = c->mtid;
    //print_read_debug(read);
    return 0;
}

void free_readmemory(struct alignedread* read) {
    free(read->readid);
    free(read->sequence);
    free(read->quality);
    free(read->barcode);
    if (read->cigs > 0) free(read->cigarlist);
}
