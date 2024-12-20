#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <htslib/thread_pool.h>
#include <iostream>
#include <cmath>
#include <math.h>

#include <pthread.h>

#include "mrand.h"
#include "Briggs.h"
#include "NtSubModels.h"
#include "RandSampling.h"
#include "getFragmentLength.h"
#include "Sampling.h"
#include "sample_qscores.h"
#include "NGSNGS_cli.h"
#include "fasta_sampler.h"
#include "add_variants.h"
#include "add_indels.h"
#include "NGSNGS_misc.h"

#define LENS 10000

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,fasta_sampler *fs,char CommandArray[LENS],const char* version){
  // Creates a header for the bamfile. The header is initialized before the function is called //
  
  char genome_len_buf[1024];
  sam_hdr_add_line(header, "HD", "VN",version, "SO", "unsorted", NULL);
  for(int i=0;i<fs->nref;i++){
    snprintf(genome_len_buf,1024,"%d", fs->seqs_l[i]);
    
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", fs->seqs_names[i], "LN", genome_len_buf, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
   
    memset(genome_len_buf,0, sizeof(genome_len_buf));
  
  }
  // Adding PG tag specifying the command used for simulations
  sam_hdr_add_pg(header,"NGSNGS","VN",version,"CL",CommandArray,NULL);
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile_nam); //outfile
}

void* ThreadInitialization(const char* version,char CommandArray[LENS],int thread_no,int threadwriteno,size_t BufferLength,
                        const char* refSseq,const char* Specific_Chr,int seed,int RandMacro,
                        const char* OutputName,outputformat_e OutputFormat,int Align,seqtype_e SeqType,int simmode,size_t reads,size_t flankingregion, const char* BedFile, int MaskBed,
                        const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,int Lowerlimit,
                        int AddAdapt,const char* Adapter_1,const char* Adapter_2,const char* Polynt,
                        int DoSeqErr,const char* QualStringFlag,int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int FixedQual,int readcycle,int readcycle_fix,
                        int doMisMatchErr,const char* SubProfile,int MisLength,const char* MisMatchMatrix,const char* M3outname,
                        float BriggsParam[4],int DoNonBiotin,int DoBiotin,int Duplicates,
                        double mutationrate, size_t referencevariations, int generations,char* VariationfileDump,
                        const char *VariantFile,int HeaderIndivIdx,const char* NameIndiv,const char* VCFfileDump,int CaptureVCF,int linkage,
                        float IndelFuncParam[4],int DoIndel,const char* IndelDumpFile, size_t genome_size, int fileAppend){


  /*
  NGSNGS overall -                        @param const char* version,char CommandArray[LENS],int thread_no,int threadwriteno,size_t BufferLength
  Reference specific -                    @param const char* refSseq,const char* Specific_Chr,int seed,int RandMacro,
  simulation and output specific -        @param const char* OutputName,outputformat_e OutputFormat,int Align,seqtype_e SeqType,int simmode,size_t reads,size_t flankingregion, const char* BedFile, int MaskBed,
  fragment length specific -              @param const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,int Lowerlimit,
  
  Additional nucleotide post 
    simulation specific -                 @param int AddAdapt,const char* Adapter_1,const char* Adapter_2,const char* Polynt,
  Sequencing error (fastq,sam) specific - @param int DoSeqErr,const char* QualStringFlag,int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int FixedQual,int readcycle,int readcycle_fix,
  Nucleotide misincorporation specific -  @param int doMisMatchErr,const char* SubProfile,int MisLength,const char* MisMatchMatrix,const char* M3outname,
  PMD specific -                          @param float BriggsParam[4],int DoNonBiotin,int DoBiotin,int Duplicates,
  Reference specific 
    stochastic variation -                @param double mutationrate, size_t referencevariations, int generations,
  Allele specific variations -            @param const char *VariantFile,int HeaderIndivIdx,const char* NameIndiv,const char* VCFfileDump,int CaptureVCF,int linkage,
  sequencing read specific stochastic 
    indels variations -                   @param float IndelFuncParam[4],int DoIndel,const char* IndelDumpFile
  */

  //creating an array with the arguments to create multiple threads;

  int nthreads=thread_no;
  pthread_t *mythreads = new pthread_t[nthreads];
  
  //allocate for reference file
  time_t t_ref = time(NULL);
  fasta_sampler *reffasta;

  int bedfilesample = 0;
  int vcfcapture = 0;
  int LD = 0;
  
  /*

    CREATION OF CONTIG TO SIMULATE SEQUENCING READS DEPENDING ON NGSNGS INPUTS

  */ 

  // Allocation for reference genomes depending on simulation type
  if(Specific_Chr == NULL && BedFile == NULL && MaskBed == 0 && CaptureVCF == 0 && linkage == 0){
    // all chromosomes/contigs/scaffolds from input reference genome
    reffasta = fasta_sampler_alloc_full(refSseq);
    fprintf(stderr,"\t-> Allocated memory for full genome with %d chromosomes/contigs/scaffolds from input reference genome (-i) with the full length %zu nt\n",reffasta->nref,reffasta->seq_l_total);
  }
  else if(Specific_Chr != NULL && BedFile == NULL && MaskBed == 0 && CaptureVCF == 0 && linkage == 0){
    // Subset of chromosomes/contigs/scaffolds from input reference genome
    reffasta = fasta_sampler_alloc_subset(refSseq,Specific_Chr);
    fprintf(stderr,"\t-> Allocated memory for subset (-chr) of the input reference genome (-i) with %d chromosomes/contigs/scaffolds with the full length %zu nt\n",reffasta->nref,reffasta->seq_l_total);
  }
  else if(Specific_Chr == NULL && BedFile != NULL && MaskBed == 0 && CaptureVCF == 0 && linkage == 0){
    // Creating pseudo chromosomes/contigs/scaffolds from input reference genome from regions of interest in bed file

    reffasta = fasta_sampler_alloc_bedentry(refSseq,BedFile,flankingregion);
    bedfilesample = 1; // to parse with threads to enable read id splitting chromosome name to create accurate coordinates
    
    fprintf(stderr,"\t-> Allocated memory for regions of interest (-incl) from the input reference genome (-i) with %d chromosomes/contigs/scaffolds with the full length %zu nt\n",reffasta->nref,reffasta->seq_l_total);
  }
  else if(Specific_Chr == NULL && BedFile != NULL && MaskBed == 1 && CaptureVCF == 0 && linkage == 0){
    // Creating pseudo chromosomes/contigs/scaffolds from input reference genome and masking specified genomic regions in bed file

    reffasta = fasta_sampler_alloc_maskbedentry(refSseq,BedFile,flankingregion);
    bedfilesample = 1;

    fprintf(stderr,"\t-> Allocated memory masking regions (-excl) within the genome from the input reference genome (-i) with a total of %d chromosomes/contigs/scaffolds with the full length %zu nt\n",reffasta->nref,reffasta->seq_l_total);
  }
  else if(Specific_Chr == NULL && BedFile == NULL && MaskBed == 0 && CaptureVCF == 1 && linkage == 0){
    // Creating pseudo chromosomes/contigs/scaffolds from those defined in both reference genome and vcf file, creating genome of interest with variations 

    if(VariantFile){
      reffasta = fasta_sampler_alloc_vcf(refSseq,VariantFile,HeaderIndivIdx,NameIndiv,flankingregion);
      vcfcapture = 1;
    }
    else{
      fprintf(stderr,"For --capture simulation, please provide -vcf file or consider simulating regions of interest (-incl) or mask regions (-excl) using bed file\n");
      exit(1);
    }
    fprintf(stderr,"\t-> Allocated memory capturing regions (--capture) with variants (-vcf) within the genome from the input reference genome (-i) with %d chromosomes/contigs/scaffolds with the full length %zu nt\n",reffasta->nref,reffasta->seq_l_total);
  }
  else if(Specific_Chr == NULL && BedFile == NULL && MaskBed == 0 && CaptureVCF == 0 && linkage == 1){
    // Creating pseudo chromosomes/contigs/scaffolds from those defined in both reference genome and vcf file, creating genome of interest with variations in LD 

    if(VariantFile){
      reffasta = fasta_sampler_alloc_vcf_LD(refSseq,VariantFile,HeaderIndivIdx,NameIndiv,flankingregion);
      LD = 1;
    }
    else{
      fprintf(stderr,"For linkage disequilibrium simulation, please provide -vcf file or consider simulating regions of interest (-incl) or mask regions (-excl) using bed file\n");
      exit(1);
    }

  }
  else{
    fprintf(stderr,"\t-> error allocation memory for %d chromosomes/contigs/scaffolds from input reference genome - please check with helppage in terms of parameters and the simulation mode\n",reffasta->nref);
    fprintf(stderr,"\t\t conflicts might arise depending on chosen -i, -chr, -bed, -vcf & -capture\n");
    exit(1);
  }
  
  fprintf(stderr, "\t-> Done reading in the reference file, walltime used =  %.2f sec\n", (float)(time(NULL) - t_ref));
  
  // rejust reads according reffasta
  fprintf(stderr, "current ref fasta genome size: %ld\n", reffasta->seq_l_total);
  reads = reads * reffasta->seq_l_total / genome_size;
  fprintf(stderr, "Adjusted Reads: %ld\n", reads);

  /*

    ALTER THE SAMPLED REFERENCE GENOME / PSEUDO CONTIGS - EITHER BIOLOGICAL OR STOCHASTIC

  */ 

  // biological variants for a given individual within a population
  if(VariantFile && CaptureVCF == 0 && linkage == 0){
    add_vcf_variants(reffasta,VariantFile,HeaderIndivIdx,NameIndiv);
    if(VCFfileDump!=NULL){
      char dumpfile1[512];
      const char* dumpfile1prefix = VCFfileDump;
      const char* dumpfile1suffix = ".fa";
      strcpy(dumpfile1,dumpfile1prefix);
      strcat(dumpfile1,dumpfile1suffix);
      const char* dumpfilefull = dumpfile1;
      dump_internal(reffasta,dumpfilefull);
    }
    fprintf(stderr, "\t-> Done adding variants from variant calling format, walltime used =  %.2f sec\n", (float)(time(NULL) - t_ref));
  }

  // stochastic varaitions representing a mutation rate
  if(mutationrate > 0.0 || referencevariations > 0){
    time_t t_mutation = time(NULL);
    mrand_t *mutation_rand = mrand_alloc(RandMacro,seed);
    const char *bases = "ACGTN";

    size_t num_variations = 0;

    if (mutationrate > 0.0){
      num_variations = (size_t) reffasta->seq_l_total*mutationrate*generations;
    }
    else{
      num_variations = referencevariations;
    }

    kstring_t *stochasticvariant;
    BGZF *bgzf_fp;
    if(VariationfileDump!=NULL){
    
      stochasticvariant =(kstring_t*) calloc(1,sizeof(kstring_t));
      stochasticvariant->s = NULL;
      stochasticvariant->l = stochasticvariant->m = 0;
      
    
      char vardumpfile[512];
      const char* vardumpfileprefix = VariationfileDump;
      const char* vardumpfilesuffix = ".txt";
      strcpy(vardumpfile,vardumpfileprefix);
      strcat(vardumpfile,vardumpfilesuffix);
      const char* vardumpfilefull = vardumpfile;

      const char* modefp2 = "wu";
      bgzf_fp = bgzf_open(vardumpfilefull,modefp2); //w
      if (!bgzf_fp) {
        fprintf(stderr, "Failed to open BGZF file: %s\n", vardumpfilefull);
        exit(1);
      }
      bgzf_mt(bgzf_fp,2,256);
    }

    for (size_t i = 0; i < num_variations;){
      int chr_idx = 0; //(int)(mrand_pop_long(mr) % (reffasta->nref));
      //Choose random chromosome index each time
      if(reffasta->nref>1)
        chr_idx = ransampl_draw2(reffasta->ws,mrand_pop(mutation_rand),mrand_pop(mutation_rand));

      long rand_val = mrand_pop_long(mutation_rand);
      size_t pos = (size_t)(abs(rand_val) % reffasta->seqs_l[chr_idx]);

      // Alter nucleotide within reference
      if (reffasta->seqs[chr_idx][pos] != 'N'){  
        char previous = reffasta->seqs[chr_idx][pos];
        char altered = bases[(int)(mrand_pop(mutation_rand)*4)]; //bases[(int)(mrand_pop_long(mr) %4)];
      
        while(previous == altered){
          altered = bases[(int)(mrand_pop(mutation_rand)*4)]; //bases[(int)(mrand_pop_long(mr) %4)];
        }
        reffasta->seqs[chr_idx][pos] = altered;
        i++;
        if(VariationfileDump!=NULL){
          ksprintf(stochasticvariant,"%s\t%lu\t%c\t%c\n",reffasta->seqs_names[chr_idx],pos,previous,altered);
        }
      }
      else{
        continue;
      }
    }

    if(VariationfileDump!=NULL){
      if (stochasticvariant->l > 0){
        assert(bgzf_write(bgzf_fp,stochasticvariant->s,stochasticvariant->l)!=0);
        stochasticvariant->l = 0;
      }

      free(stochasticvariant->s);
      free(stochasticvariant);

      // close the output files
      if(bgzf_fp!=NULL){
        bgzf_close(bgzf_fp);
      }    
    }

    free(mutation_rand);
    fprintf(stderr, "\t-> Done adding %zu stochastic variants to reference genome, walltime used =  %.2f sec\n", num_variations,(float)(time(NULL) - t_mutation));
  }

  /*

    FOLLOWING THE ALLOCATION,CREATION AND ALTERATION OF THE REFERENCE GENOME - CONTINUE WITH CREATOMG POTENTIAL OUTPUT FILES AND DISTRIBUTIONS USED FOR SAMPLING

  */ 

  if (reffasta->seqs != NULL){
    Parsarg_for_Sampling_thread *struct_for_threads = new Parsarg_for_Sampling_thread[nthreads];

    // declare files and headers
    BGZF **bgzf_fp = (BGZF **) calloc(3,sizeof(BGZF *));

    samFile *SAMout = NULL;
    sam_hdr_t *SAMHeader = NULL;
    htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
    htsThreadPool p = {NULL, 0};

    char file1[512];
    char file2[512];
    const char* fileprefix = OutputName;
    strcpy(file1,fileprefix);
    strcpy(file2,fileprefix);

    const char* suffix1 = NULL;
    const char* suffix2 = NULL;
    const char* mode = NULL;
    int alnformatflag = 0;

    switch(OutputFormat){
      // determine output prefix and for output format and simulation mode
      case faT:
        mode = "wu";
        if(SE==SeqType)
          suffix1 = ".fa";
        else{
          suffix1 = "_R1.fa";
          suffix2 = "_R2.fa";
        }
        break;
      case fagzT:
        mode = "wb";
        if(SE== SeqType)
          suffix1 = ".fa.gz";
        else{
          suffix1 = "_R1.fa.gz";
          suffix2 = "_R2.fa.gz";
        }
        break;
      case fqT:
        mode = "wu";
        if(SE ==SeqType)
          suffix1 = ".fq";
        else{
          suffix1 = "_R1.fq";
          suffix2 = "_R2.fq";
        }
        break;
      case fqgzT:
        if(fileAppend==1) mode = "a";
        else mode = "w";
        if(SE==SeqType)
          suffix1 = ".fq.gz";
        else{
          suffix1 = "_R1.fq.gz";
          suffix2 = "_R2.fq.gz";
        }
        break;

      case samT:
        mode = "ws";
        suffix1 = ".sam";
        alnformatflag++;
        break;

      case bamT:
        mode = "wb";
        suffix1 = ".bam";
        alnformatflag++;
        break;
      case cramT:
        mode = "wc";
        suffix1 = ".cram";
        alnformatflag++;
        break;
      default:
        fprintf(stderr,"\t-> Fileformat is currently not supported \n");
        break;
    }
     
    strcat(file1,suffix1);

    fprintf(stderr,"\t-> File output name is %s\n",file1);
    const char* filename1 = file1;
    const char* filename2 = NULL;

    if(alnformatflag == 0){
      // for fasta and fastq formats simply store the files
      int mt_cores = threadwriteno;
      int bgzf_buf = 256;
      
      bgzf_fp[0] = bgzf_open(filename1,mode);
      bgzf_mt(bgzf_fp[0],mt_cores,bgzf_buf);
      
      if(PE==SeqType){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp[1] = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp[1],mt_cores,bgzf_buf);
      }
    }
    else{
      // Create sam file header 
      char *ref =(char*) malloc(strlen(".fasta.gz") + strlen(refSseq) + 2);
      sprintf(ref, "reference=%s", refSseq);
      
      // Save reference file name for header creation of the sam output
      // hts_opt_add((hts_opt **)&fmt_hts->specific,ref);
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();

      if(threadwriteno>0){
        if (!(p.pool = hts_tpool_init(threadwriteno))) {
          fprintf(stderr, "Error creating thread pool\n");
          exit(1);
        }
        hts_set_opt(SAMout, HTS_OPT_THREAD_POOL, &p);
      }
      hts_set_opt(SAMout, CRAM_OPT_REFERENCE, refSseq);
      // generate header
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,reffasta,CommandArray,version);

      free(ref);
      // hts_opt_free((hts_opt *)fmt_hts->specific);
    }

    //Read in the fragment length file (-lf) before creating threads

    int no_elem;double* Frag_freq;int* Frag_len;
    if(SizeDistType==1){
        Frag_len = new int[LENS];Frag_freq = new double[LENS];
        ReadLengthFile(no_elem,Frag_len,Frag_freq,Sizefile);
    }
    else{no_elem = -1;}

    // prepare information for quality specific simulation information
    int inferred_readcycle = 0;
    const char *freqfile_r1;
    const char *freqfile_r2;
    int outputoffset = qualstringoffset;
    ransampl_ws ***QualDist = NULL;
    char nt_qual_r1[1024];
    ransampl_ws ***QualDist2 = NULL;
    char nt_qual_r2[1024];
    double ErrArray_r1[1024];
    double ErrArray_r2[1024];
    freqfile_r1 = QualProfile1;

    // create sequence read quality distributions
    if(strcasecmp("true",QualStringFlag)==0){
      if(QualProfile1 != NULL && FixedQual == 0){
        QualDist = ReadQuality(nt_qual_r1,ErrArray_r1,outputoffset,freqfile_r1,inferred_readcycle);
        if(PE==SeqType){
          freqfile_r2 = QualProfile2;
          QualDist2 = ReadQuality(nt_qual_r2,ErrArray_r2,outputoffset,freqfile_r2,inferred_readcycle);
        }
      }
    }
    
    int maxsize = 20;
    char polynucleotide;

    if (Polynt != NULL && strlen(Polynt) == 1){polynucleotide = (char) Polynt[0];}
    else{polynucleotide = 'F';}

    //generating mismatch matrix to parse for each string either from -mf or -m3 
    double* MisMatchFreqArray = new double[LENS];
    int mismatchcyclelength = 0;
    int numElements = 0;
    if (SubProfile != NULL){
      // nucleotide subsitution file
      MisMatchFileArray(MisMatchFreqArray,SubProfile,mismatchcyclelength,numElements);
      //exit(1);
    }
    if (MisMatchMatrix != NULL){
      // bdamage.gz file from metaDMG-cpp (05-09-2024)
      const char* dumpM3full = NULL;
      if(M3outname!=NULL){
        char dumpM3[512];
        const char* dumpM3prefix = M3outname;
        const char* dumpM3suffix = ".txt";
        strcpy(dumpM3,dumpM3prefix);
        strcat(dumpM3,dumpM3suffix);
        dumpM3full = dumpM3;
      }
      MisMatchMetaFileArray(MisMatchFreqArray,MisMatchMatrix,mismatchcyclelength,numElements,dumpM3full);  
    }
    // output file for internal information from stochastic indels
    if(IndelDumpFile!=NULL){
      char IndelFile[512];
      const char* IndelSuffix = ".txt";
      const char* IndelPrefix = IndelDumpFile;
      strcpy(IndelFile,IndelPrefix);
      strcat(IndelFile,IndelSuffix);
      const char* modefp2 = "wu";
      bgzf_fp[2] = bgzf_open(IndelFile,modefp2); //w
      bgzf_mt(bgzf_fp[2],threadwriteno,256); //
    }

    /*

      CREATE SAMPLING THREADS WITH ALL INFORMATION REQUIRED FOR SIMULATING SEQUENCES

    */ 

    for (int i = 0; i < nthreads; i++){
      struct_for_threads[i].reffasta = reffasta;
  
      // The output format, output files, and structural elements for SAM outputs
      struct_for_threads[i].OutputFormat = OutputFormat;
      struct_for_threads[i].SeqType = SeqType;
      struct_for_threads[i].bgzf_fp = bgzf_fp;
      struct_for_threads[i].SAMout = SAMout;
      struct_for_threads[i].SAMHeader = SAMHeader;
      struct_for_threads[i].LengthData = 0;
      struct_for_threads[i].MaximumLength = maxsize;
      struct_for_threads[i].list_of_reads = (bam1_t**) malloc(sizeof(bam1_t)*maxsize); // need to free this space
      for(int j=0; j<maxsize;j++){struct_for_threads[i].list_of_reads[j]=bam_init1();} // but also destroy the bam_init1 objects    

      // Thread generation and sampling specific information
      struct_for_threads[i].threadno = i;
      struct_for_threads[i].totalThreads = nthreads;
      struct_for_threads[i].threadseed = seed;
      struct_for_threads[i].rng_type = RandMacro;
      struct_for_threads[i].simmode = simmode;
      struct_for_threads[i].bedfilesample = bedfilesample;
      struct_for_threads[i].VCFcapture = vcfcapture;

      // Sequence alteration models
        // 1) nucleotide quality score and sequencing errors,  
      struct_for_threads[i].QualFlag = QualStringFlag;
      struct_for_threads[i].DoSeqErr = DoSeqErr;
      struct_for_threads[i].NtQual_r1 = nt_qual_r1;
      struct_for_threads[i].NtQual_r2 = nt_qual_r2;
      struct_for_threads[i].QualDist_r1 = QualDist;
      struct_for_threads[i].QualDist_r2 = QualDist2;
      struct_for_threads[i].FixedQual_r1r2 = FixedQual;
      struct_for_threads[i].qsreadcycle = (int) readcycle_fix;

        // 2) Stochastic indel
      struct_for_threads[i].NtErr_r1 = ErrArray_r1;
      struct_for_threads[i].NtErr_r2 = ErrArray_r2;
      struct_for_threads[i].maxreadlength = (int) readcycle;
      struct_for_threads[i].IndelFuncParam = IndelFuncParam;
      struct_for_threads[i].DoIndel = DoIndel;
      struct_for_threads[i].IndelDumpFile = IndelDumpFile;
      
        // 3) PMD models
      struct_for_threads[i].DoNonBiotin = DoNonBiotin;
      struct_for_threads[i].DoBiotin = DoBiotin;
      struct_for_threads[i].BriggsParam = BriggsParam;
      struct_for_threads[i].Duplicates = Duplicates;

        // 4) mismatch matrices either from -mf or -m3 with same internal format
      struct_for_threads[i].MisMatch = MisMatchFreqArray;
      struct_for_threads[i].doMisMatchErr = doMisMatchErr;
      struct_for_threads[i].MisLength = (int) mismatchcyclelength;

      // Fragment lengths 
      struct_for_threads[i].FragLen = Frag_len;
      struct_for_threads[i].FragFreq = Frag_freq;
      struct_for_threads[i].No_Len_Val = no_elem;
      struct_for_threads[i].FixedSize = FixedSize;
      struct_for_threads[i].distparam1 = val1;
      struct_for_threads[i].distparam2 = val2;
      struct_for_threads[i].LengthType = SizeDistType;
      struct_for_threads[i].lowerlimit = Lowerlimit;

      // Sequence output specific
      struct_for_threads[i].BufferLength = BufferLength;

      // Additional information for sequence reads
      struct_for_threads[i].AddAdapt = AddAdapt;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      struct_for_threads[i].Adapter_2 = Adapter_2;
      struct_for_threads[i].PolyNt = polynucleotide;
      struct_for_threads[i].Align = Align;

    }

    size_t ThreadReads = (size_t) floor( reads / (double) thread_no);
    for (int i = 0; i < nthreads-1; i++){
      struct_for_threads[i].reads = ThreadReads;
    }
    struct_for_threads[nthreads-1].reads = reads - (ThreadReads*(nthreads-1));

    // initialize and create the threads
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    if(nthreads==1){
      Sampling_threads(struct_for_threads);
    }
    else{
      for (int i = 0; i < nthreads; i++){
	      int ret = pthread_create(&mythreads[i],&attr,Sampling_threads,&struct_for_threads[i]);
        if (ret != 0) {
          fprintf(stderr, "Error creating thread: %s\n");
          exit(1);
        }
      }

      for (int i = 0; i < nthreads; i++){  
	      pthread_join(mythreads[i],NULL);
      }
    }
    
    // close the output files
    if(bgzf_fp[0]!=NULL){
      bgzf_close(bgzf_fp[0]);
    }
    if(bgzf_fp[1]!=NULL){
      bgzf_close(bgzf_fp[1]);
    }
    if(bgzf_fp[2]!=NULL){
      bgzf_close(bgzf_fp[2]);
    }
    free(bgzf_fp); //free the calloc
     
    if(SAMHeader)
    sam_hdr_destroy(SAMHeader);
    if(SAMout)
      sam_close(SAMout);
    if (p.pool)
      hts_tpool_destroy(p.pool);
    
    if(CaptureVCF == 0 && linkage == 0){
      fasta_sampler_destroy(reffasta);
    }
    else{
      fasta_sampler_destroy_captureLD(reffasta);
    }

    // clean up memory
    for(int i=0;i<nthreads;i++)
      free(struct_for_threads[i].list_of_reads);

    
    delete[] mythreads;
    if(QualProfile1 != NULL && FixedQual == 0){
      for(int base=0;base<5;base++){
        for(int pos = 0 ; pos< (int) readcycle;pos++){
          ransampl_free(QualDist[base][pos]);
        }
        delete[] QualDist[base];
      }
      delete[] QualDist;

      if(PE==SeqType){
        for(int base=0;base<5;base++){
          for(int pos = 0 ; pos< (int) readcycle;pos++){
            ransampl_free(QualDist2[base][pos]);
          }
          delete[] QualDist2[base];
        }
        delete[] QualDist2;
      }
    }
    
    free(fmt_hts);
    if(SizeDistType==1){
      delete[] Frag_freq;
      delete[] Frag_len;
    }

    delete[] struct_for_threads;

    delete[] MisMatchFreqArray;
    
    fflush(stderr);
  }
  return NULL;
}

