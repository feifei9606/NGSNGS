#include "../mrand.h"
#include "../fasta_sampler.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <iostream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>

typedef struct{
  const char *bcffilename;
  const char *bedfilename;
  int range;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Stochastic variation on input reference genome\n");
  fprintf(fp,"Usage\n./VCFtoBED -vcf <vcf_file> -r <range before and after vcf position> -o <Output file name for BED>\n");
  fprintf(fp,"\nExample\n./VCFtoBED -vcf ../Test_Examples/ChrMtSubDeletionDiploid.vcf -bed lol.bed -r 10\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-vcf | --vcfin: \t\t The input vcf file containing positions with variants\n");
  fprintf(fp,"-bed | --bedout: \t\t The output name for the created bed file with only regions present in the vcf file\n");
  fprintf(fp,"-r | --range: \t\t Integer value specifying the number of nucleotides before and after the vcf position, default = 10\n");
  exit(1);
  return 0;
}
 

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->bcffilename = NULL;
  mypars->bedfilename = NULL;
  mypars->range = 10;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-vcf",*argv)==0 || strcasecmp("--vcfin",*argv)==0){
      mypars->bcffilename = strdup(*(++argv));
    }
    else if(strcasecmp("-bed",*argv)==0 || strcasecmp("--bedout",*argv)==0){
      mypars->bedfilename = strdup(*(++argv));
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--range",*argv)==0){
      mypars->range = atoi(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}

int VCFtoBED(int argc,char **argv){
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                          strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else{
    mypars = getpars(argc,argv);
    const char *bcffilename = mypars->bcffilename;
    const char *bedfilename = mypars->bedfilename;
    int range = mypars->range;
    
    htsFile *bcf = bcf_open(bcffilename, "r");
    if (!bcf) {
      fprintf(stderr, "Error: Unable to open BCF file: %s\n", bcffilename);
      return 1;
    }

    bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
    bcf1_t *brec = bcf_init();

    FILE *bedfile = fopen(bedfilename, "w");
    if (!bedfile) {
      fprintf(stderr, "Error: Unable to open BED file for writing: %s\n", bedfilename);
      return 1;
    }

    int seqnames_l;
    const char **bcf_chrom_names = bcf_hdr_seqnames(bcf_head,&seqnames_l);
    
    bcf_idpair_t *chr_info = bcf_head->id[BCF_DT_CTG];
    
    //for (int i = 0; i < bcf_head->n[BCF_DT_CTG]; ++i)
    //  fprintf(stderr,"%s\t%d\n", ctg[i].key, ctg[i].val->info[0]);


    while (bcf_read(bcf, bcf_head, brec) == 0) {
      bcf_unpack((bcf1_t *)brec, BCF_UN_ALL);
      int fai_chr = brec->rid;
      if (fai_chr == -1) {
        fprintf(stderr, "Error: Chromosome name not found in the reference\n");
        return 1;
      }
      int chr_end = chr_info[fai_chr].val->info[0];
      if((brec->pos - range) < 1){
        fprintf(bedfile, "%s\t%lld\t%lld\n", bcf_hdr_id2name(bcf_head, brec->rid), 1, brec->pos + range);
      }
      else if((brec->pos + range) > chr_end){
        fprintf(bedfile, "%s\t%lld\t%lld\n", bcf_hdr_id2name(bcf_head, brec->rid), brec->pos - range, chr_end);
      }
      else{
        fprintf(bedfile, "%s\t%lld\t%lld\n", bcf_hdr_id2name(bcf_head, brec->rid), brec->pos - range, brec->pos + range);
      }
    }

    fclose(bedfile);
    bcf_hdr_destroy(bcf_head);
    bcf_destroy(brec);
    bcf_close(bcf);
  }

  return 0;
}


#ifdef __WITH_MAIN__
int main(int argc,char **argv){
  VCFtoBED(argc,argv);
  return 0;
}
#endif

/*
g++ VCFtoBED.cpp -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto /willerslev/users-shared/science-snm-willerslev-wql443/WP1/htslib/libhts.a -D __WITH_MAIN__ -o ChrRandVar

./VCFtoBED -vcf ../Test_Examples/ChrMtSubDeletionDiploid.vcf -bed lol.bed -r 10
*/