#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"
#include "Sampling.h"

typedef struct{
  int threads;
  size_t reads;
  int Glob_seed;
  const char *OutFormat;
  const char *OutName;
  const char *Reference;
  const char *Seq;
  const char *Adapter1;
  const char *Adapter2;
  const char *QualProfile1;
  const char *QualProfile2;
  const char *Briggs;
  int Length;
  const char *LengthFile;
}argStruct;

float myatof(char *str){
  if (str==NULL)
    fprintf(stderr,"Could not parse Briggs parameters, provide <nv,Lambda,Delta_s,Delta_d>");
  
  return atof(str);
}

void Sizebreak(char *str){
  fprintf(stderr,"Could not parse the length parameters, provide either fixed length size (-l) \n or parse length distribution file (-lf)");
}

int HelpPage(FILE *fp){
  fprintf(fp,"Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 \n\n");
  fprintf(fp,"Default usage: ./ngsngs <input reference> <numer of reads> \t Generates single-end fasta file named output\n\n");
  fprintf(fp,"Usage: ./ngsngs <options> <input reference> <numer of reads> <output file> \n\n");
  fprintf(fp,"Required options: \n");
  fprintf(fp,"-i | --input: \t\t Reference file in .fasta format to sample reads from\n");
  fprintf(fp,"-r | --reads: \t\t Number of reads to simulate\n");
  fprintf(fp,"\n");
  fprintf(fp,"Optional options: \n");
  fprintf(fp,"-h | --help: \t\t Print help page\n");
  fprintf(fp,"-v | --version: \t\t Print help page\n");
  fprintf(fp,"-o | --output: \t\t Prefix of output file name, with default extension in fasta format (.fa)\n");
  fprintf(fp,"-f | --format: \t\t File format of the simulated outpur reads\n");
  fprintf(fp,"\t <.fa||.fasta>\t nucletide sequence \n \t <.fa.gz||.fasta.gz>\t compressed nucletide sequence \n \t <.fq||.fastq>\t nucletide sequence with corresponding quality score \n \t <.fq.gz||.fastq.gz>\t compressed nucletide sequence with corresponding quality score \n \t <.sam||bam>\t Sequence Alignment Map format\n");
  fprintf(fp,"-t | --threads: \t Number of threads to use for simulation\n");
  fprintf(fp,"-s | --seed: \t\t Random seed, default value being computer time\n");
  fprintf(fp,"-a1 | --adapter1: \t Adapter sequence to add for simulated reads (SE) or first read pair (PE)\n");
  fprintf(fp,"-a2 | --adapter2: \t Adapter sequence to add for second read pair (PE) \n");
  fprintf(fp,"-seq | --sequencing: \t Simulate single-end or paired-end\n");
  fprintf(fp,"\t <SE>\t single-end \n \t <PE>\t paired-end\n");
  fprintf(fp,"-b | --briggs: <nv,Lambda,Delta_s,Delta_d>\t Parameters for the damage patterns (Briggs et al., 2007)\n");
  fprintf(fp,"\t <nv>\t Nick rate pr site \n \t <Lambda>\t Geometric distribution parameter for overhang length\n \t <Delta_s>\t PMD rate in single-strand regions\n \t <Delta_s>\t PMD rate in double-strand regions");
  exit(1);
  return 0;
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->threads = 1;
  mypars->reads = -1;
  mypars->Glob_seed = (int) time(NULL);
  mypars->OutFormat = "fa";
  mypars->OutName = "output";
  mypars->Seq = "SE";
  mypars->Reference = NULL;
  mypars->Adapter1 = NULL;
  mypars->Adapter2 = NULL;
  mypars->QualProfile1 = NULL;
  mypars->QualProfile2 = NULL;
  mypars->Briggs = NULL; //"0.024,0.36,0.68,0.0097";
  mypars->Length = -1;
  mypars->LengthFile = NULL;
  while(*argv){
    // Required
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Reference = strdup(*(++argv));
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      mypars->reads = atoi(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->OutName = strdup(*(++argv));
    }
    // Optional
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->threads = atoi(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Glob_seed = atoi(*(++argv));
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      mypars->Seq = strdup(*(++argv));
      if(strcasecmp("SE",mypars->Seq)!=0 && strcasecmp("PE",mypars->Seq)!=0){HelpPage(stderr);} 
    }
    else if(strcasecmp("-a1",*argv)==0 || strcasecmp("--adapter1",*argv)==0){
      mypars->Adapter1 = strdup(*(++argv));
    }
    else if(strcasecmp("-a2",*argv)==0 || strcasecmp("--adapter2",*argv)==0){
      mypars->Adapter2 = strdup(*(++argv));
    }
    else if(strcasecmp("-q1",*argv)==0 || strcasecmp("--quality1",*argv)==0){
      mypars->QualProfile1 = strdup(*(++argv));
    }
    else if(strcasecmp("-q2",*argv)==0 || strcasecmp("--quality2",*argv)==0){
      mypars->QualProfile2 = strdup(*(++argv));
    }
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      mypars->OutFormat = strdup(*(++argv));
    }
    else if(strcasecmp("-b",*argv)==0 || strcasecmp("--briggs",*argv)==0){
      mypars->Briggs = strdup(*(++argv)); //double nv, double lambda, double delta_s, double delta
      // "0.01,0.036,1.02,2.0"
    }
    else if(strcasecmp("-l",*argv)==0 || strcasecmp("--length",*argv)==0){
      mypars->Length = atoi(*(++argv));
    }
    else if(strcasecmp("-lf",*argv)==0 || strcasecmp("--lengthfile",*argv)==0){
      mypars->LengthFile = strdup(*(++argv));
    }
    
    /*else{
      fprintf(stderr,"unrecognized input option, see NGSNGS help page\n\n");
      HelpPage(stderr);
      } */
    
    // -e1 +2 || --error1 +2 
    ++argv;
  }
  return mypars;
}

// ------------------------------ //

int main(int argc,char **argv){
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                        strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else{
    mypars = getpars(argc,argv);
    clock_t t = clock();
    time_t t2 = time(NULL);

    const char *fastafile = mypars->Reference; //"/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa";
    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(fastafile);

    fprintf(stderr,"\t-> fasta load \n");
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);
    int Glob_seed = mypars->Glob_seed; //(int) time(NULL);
    int threads = mypars->threads;
    size_t No_reads = mypars->reads; //1e1;

    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
    fprintf(stderr,"\t-> Seed used: %d\n",Glob_seed);
    fprintf(stderr,"\t-> Number of threads used: %d \n",threads);
    fprintf(stderr,"\t-> Number of simulated reads: %zd\n",No_reads);

    const char* Adapt_flag;
    const char* Adapter_1;
    const char* Adapter_2;
    if (mypars->Adapter1 != NULL){
      Adapt_flag = "true";
      Adapter_1 = mypars->Adapter1;
      Adapter_2 = mypars->Adapter2;
      //Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
    }
    else{Adapt_flag = "false";}

    const char* OutputFormat = mypars->OutFormat;
    
    const char* QualProfile1; const char* QualProfile2;
    QualProfile1 = mypars->QualProfile1; QualProfile2 = mypars->QualProfile2;
    if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0 || strcasecmp("bam",OutputFormat)==0){
      if (QualProfile1 == NULL){
        fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
        exit(0);
      }
      if(strcasecmp("PE",mypars->Seq)==0 && mypars->QualProfile2 == NULL){
        fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
        exit(0);
      }
    }
    
    int qualstringoffset = 0;
    if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){qualstringoffset = 33;}

    int Thread_specific_Read = static_cast<int>(No_reads/threads);

    const char* filename = mypars->OutName; //"chr22_out";
    const char* Seq_Type = mypars->Seq;
    
    int FixedSize = mypars->Length;
    const char* Sizefile = mypars->LengthFile;

    if (Sizefile == NULL){
      if (FixedSize == -1){
        fprintf(stderr,"Could not parse the length parameters, provide either fixed length size (-l) or parse length distribution file (-lf) see helppage (-h).");
        exit(0);
      }
    }
    
    const char* Briggs_Flag;
    float Param[4];
    if (mypars->Briggs != NULL){
      char* BriggsParam = strdup(mypars->Briggs);
      Param[0] = myatof(strtok(BriggsParam,"\", \t"));
      Param[1] = myatof(strtok(NULL,"\", \t"));
      Param[2] = myatof(strtok(NULL,"\", \t"));
      Param[3] = myatof(strtok(NULL,"\", \t"));
      Briggs_Flag = "True";
    }
    else{Briggs_Flag = "False";}

    Create_se_threads(seq_ref,threads,Glob_seed,Thread_specific_Read,filename,
                      Adapt_flag,Adapter_1,Adapter_2,OutputFormat,Seq_Type,
                      Param,Briggs_Flag,Sizefile,FixedSize,qualstringoffset,QualProfile1,QualProfile2);

    fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }

  /*delete mypars; 
   200 (96 direct, 104 indirect) bytes in 1 blocks are definitely lost in loss record 4 of 4
   operator new(unsigned long) (vg_replace_malloc.c:344)
   getpars(int, char**) (in /home/wql443/WP1/NGSNGS/a.out)

   adding the delete will give me 
   bytes in 1 blocks are definitely lost in loss record ->  strdup (in /usr/lib64/libc-2.17.so)
  */
}


// g++ NGSNGS_func.cpp atomic_fq.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// valgrind --tool=memcheck --leak-check=full  --track-origins=yes ./a.out
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -f bam -o chr22
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 1000 -f fq -s 1 -o chr22_out
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq SE -f fq -o chr22
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq PE -f fq -o chr22
//  make HTSSRC=/home/wql443/scratch/htslib/
// ./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq PE -f fq -o chr22

//g++ NGSNGS_func.cpp Sampling.cpp ArgParse.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl