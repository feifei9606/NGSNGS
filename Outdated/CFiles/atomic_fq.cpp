#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"

char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],int chr_size_cumm[]){
  
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){

    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    //sprintf(&genome[strlen(genome)],data);
    //strcat(genome,data);  //Both gives conditional jump or move error
    if (data != NULL){
      sprintf(genome+strlen(genome),data); 
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total){
  // Creates a header for the bamfile. The header is initialized before the function is called //

  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
    
  // Creating header information
  char *name_len_char =(char*) malloc(1024);
  for(int i=0;i<chr_total;i++){
    const char *name = faidx_iseq(seq_ref,i);
    int name_len =  faidx_seq_len(seq_ref,name);
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", name_len_char, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
  }
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);
  //free(name_len_char);
}

pthread_mutex_t Fq_write_mutex;

// ---------------------- SINGLE-END ---------------------- //
struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  char *genome; // The actual concatenated genome
  int chr_no;
  int threadno;
  int *size_cumm;
  const char **names;

  int* FragLen;
  double* FragFreq;
  int No_Len_Val;

  double* Qualfreq;
  int threadseed;
  size_t reads;
  
  BGZF *bgzf;
  samFile *SAMout;
  sam_hdr_t *SAMHeader;
  bam1_t **list_of_reads;
  int l;
  int m;
  const char* Adapter_flag;
  const char* Adapter_1;

  const char* OutputFormat;

};
      
void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;
  //std::cout << "se run "<< std::endl;
  time_t t4=time(NULL);
  // creating random objects for all distributions.
  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno;
  
  size_t genome_len = strlen(struct_obj->genome);

  //coverage method2
  int fq_ascii_offset = 33;
  char seq_r1[1024] = {0};
  char seq_r1_mod[1024] = {0};
  char read[1024] = {0};
  char readadapt[1024] = {0};
  // for the coverage examples
  int reads = struct_obj -> reads;

  //float cov_current = 0;
  size_t rand_start;
  //int nread = 0;

  char qual[1024] = "";
  //int D_i = 0;
  int localread = 0;
  int iter = 0;
  int current_reads_atom = 0;
  //while (current_cov_atom < cov) {
  //while (current_reads_atom < reads){
  unsigned int test = 0;
  while (current_reads_atom < reads){
    //fprintf(stderr,"While loop \n");
    double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX);
    double rand_val2 = rand_val*RAND_MAX;
    unsigned int test = (unsigned int) rand_val2;
    double rand_val3 = myrand((unsigned int) rand_val2); //((double) rand_r(&test)/ RAND_MAX);// 
    //fprintf(stderr,"%lf \t %lf \t %u \t %lf \n",rand_val,rand_val2,test,rand_val3);
    //fprintf(stderr,"%lf \t %lf \n",rand_val,rand_val3);
    rand_start = rand_val3 * (genome_len-300); //genome_len-100000;
    //fprintf(stderr,"Start 1 %ld \t start 2 %ld \n",rand_start,(size_t) (rand_val * (genome_len-300)));
    int lengthbin = BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, rand_val);
    //fprintf(stderr,"%d\n",lengthbin);//BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, 0.083621));
    int fraglength = struct_obj->FragLen[lengthbin];//75; //struct_obj->FragLen[lengthbin];
    
    //fprintf(stderr,"random val %lf, start %d, fraglength %d\n",rand_val,rand_start,fraglength);
    //if (rand_start < 30000){rand_start += 30000;}
    
    //fprintf(stderr,"%d\n",rand_start);//struct_obj->FragLen[6]);
    //fprintf(stderr,"RANDOM START LENGTH %d \n",rand_start);
    //identify the chromosome based on the coordinates from the cummulative size array
    
    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}

    //deamination should be here!
    SimBriggsModel(seq_r1, seq_r1_mod, fraglength, 0.024, 0.36, 0.68, 0.0097,loc_seed);

    // case 1
    if (fraglength > 150){strncpy(seq_r1,struct_obj->genome+rand_start-1,150);}
    // case 2
    else if (fraglength <= 150){strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);}

    int rand_id = rand_val * fraglength-1; //100
    
    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seq_r1,'N');
    pch2 = strrchr(seq_r1,'N');
    //if (pch != NULL){continue;}
    if ((int )(pch-seq_r1+1) == 1 && (int)(pch2-seq_r1+1)  == strlen(seq_r1)){memset(seq_r1, 0, sizeof seq_r1);}
    else{
      int seqlen = strlen(seq_r1);

      //for (int j = 0; j < fraglength; j++){D_total += 1;}
      //std::time(nullptr)
      // SimBriggsModel(seq_r1, seq_r1_mod, fraglength, 0.024, 0.36, 0.68, 0.0097,loc_seed);
      
      int strand; 
      
      if (struct_obj->OutputFormat != "bam"){strand = (int) rand_r(&loc_seed)%2;}//1;//rand() % 2;

      //fprintf(stderr,"----------------------------\n strand1 Sequence r1 before %s \n ",seq_r1);
      if (strand == 0){
        DNA_complement(seq_r1);
        reverseChar(seq_r1);
      }
      //fprintf(stderr,"strand1 Sequence r1 after  %s \n ----------------------------\n",seq_r1);    
      
      char READ_ID[1024]; int read_id_length;
      read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%s:%d-%d_length:%d", struct_obj->threadno, rand_id,strand,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength);
      //fprintf(stderr,"%s\n",READ_ID);
      if (struct_obj->Adapter_flag == "true"){
        strcpy(read, seq_r1);
        strcat(read,struct_obj->Adapter_1);
        strncpy(readadapt, read, 150);
        if (struct_obj -> OutputFormat == "fa" || struct_obj -> OutputFormat == "fa.gz"){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,readadapt);
        }
        if (struct_obj -> OutputFormat == "fq" || struct_obj -> OutputFormat == "fq.gz"){
          Read_Qual_new(readadapt,qual,loc_seed,struct_obj->Qualfreq,33);
            
          ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,readadapt,qual);
        }
        if (struct_obj -> OutputFormat == "sam" || struct_obj -> OutputFormat == "bam"){
          Read_Qual_new(readadapt,qual,loc_seed,struct_obj->Qualfreq,0);
          // INSERT BAM PART HERE
        }
      }
      else{
        if (struct_obj -> OutputFormat == "fa" || struct_obj -> OutputFormat == "fa.gz"){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,seq_r1);

        }
        if (struct_obj -> OutputFormat == "fq" || struct_obj -> OutputFormat == "fq.gz"){
          Read_Qual_new(seq_r1,qual,loc_seed,struct_obj->Qualfreq,33);
          ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,seq_r1,qual);

        }
        if (struct_obj -> OutputFormat == "sam" || struct_obj -> OutputFormat == "bam"){
          Read_Qual_new(seq_r1,qual,loc_seed,struct_obj->Qualfreq,0);
          //fprintf(stderr," TRYING THE BAM OUTPUT %s \n",seq_r1);
          ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);

          /*pthread_mutex_lock(&Fq_write_mutex);
          bam_set1(bam_file_chr,read_id_length,READ_ID,0,chr_idx,0,0,n_cigar,cigar,-1,-1,0,0,seq_r1,qual,l_aux);
          sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,bam_file_chr);
          pthread_mutex_unlock(&Fq_write_mutex);*/

          // fprintf(stderr,"else statemetn \n");
          // INSERT BAM PART HERE
        }
      }

      if (struct_obj -> OutputFormat == "bam"){
        //fprintf(stderr,"Inside outputformat loop");
        size_t n_cigar;
        uint32_t cigar_bitstring; const uint32_t *cigar;
        size_t l_aux = 0; uint8_t mapq = 60;
        hts_pos_t min_beg, max_end, insert;
        cigar_bitstring = bam_cigar_gen(strlen(struct_obj->fqresult_r1->s), BAM_CMATCH); 
        n_cigar = 1; // Number of cigar operations, 1 since we only have matches
        uint32_t cigar_arr[] = {cigar_bitstring}; //converting uint32_t {aka unsigned int} to const uint32_t* 
        cigar = cigar_arr;
        min_beg = rand_start-struct_obj->size_cumm[chr_idx] - 1;
        uint16_t flag;
        if (strand == 0){flag = 16;}
        else{flag = 0;}

        bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag,chr_idx,min_beg,mapq,n_cigar,cigar,-1,-1,0,strlen(seq_r1),seq_r1,qual,l_aux);

        if (struct_obj->l < struct_obj->m){   
          pthread_mutex_lock(&Fq_write_mutex);
          for (int k = 0; k < struct_obj->l; k++){
            sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]);
          }
          //fprintf(stderr,"\n sam_write works\n");
          struct_obj->l = 0;
          pthread_mutex_unlock(&Fq_write_mutex);
        }

        struct_obj->fqresult_r1->l =0;
      }

      if (struct_obj -> OutputFormat == "fq"){
        //fprintf(stderr,"\t Buffer mutex with thread no %d\n", struct_obj->threadno);fflush(stderr);
        if (struct_obj->fqresult_r1->l > 30000000){
          pthread_mutex_lock(&Fq_write_mutex);
          bgzf_write(struct_obj->bgzf,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->fqresult_r1->l =0;
        }
      }
    
      memset(qual, 0, sizeof(qual));  
      memset(seq_r1, 0, sizeof seq_r1);
      memset(seq_r1_mod, 0, sizeof seq_r1_mod);
      chr_idx = 0;
      //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
      iter++;
      localread++;
      current_reads_atom++;
    }
  }
  if (struct_obj->fqresult_r1->l > 0){
    //fprintf(stderr,"\t last Buffer mutex with thread no %d\n", struct_obj->threadno);fflush(stderr);
    if (struct_obj -> OutputFormat == "fq"){
      pthread_mutex_lock(&Fq_write_mutex);
      bgzf_write(struct_obj->bgzf,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
      pthread_mutex_unlock(&Fq_write_mutex);
      struct_obj->fqresult_r1->l =0;
    }
  } 

  //delete[] struct_obj->Qualfreq;
  //delete[] struct_obj->sizearray;
  //consider freeing these after the join operator
  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  fprintf(stderr,"\t number of reads generated by this thread %d \n",localread);
  //fprintf(stderr, "\t[ALL done] walltime spend in thread %d =  %.2f sec\n", struct_obj->threadno, (float)(time(NULL) - t4));  
  //std::cout << "thread done" << std::endl;
  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, int reads,const char* Adapt_flag,const char* Adapter_1,const char* OutputFormat){
  time_t t3=time(NULL);
  //creating an array with the arguments to create multiple threads;
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];
  
  int chr_total = faidx_nseq(seq_ref);
  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_size_cumm[chr_total+1];
  
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);

  if (genome_data != NULL){
    fprintf(stderr,"\t-> Full genome function run!\n");
    fprintf(stderr,"\t-> Full genome size %lu \n",strlen(genome_data));
    
    //std::cout << " genome length " << genome_len << std::endl;
    Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];
 
    // creating file type and name
    char file[80];
    const char* fileprefix = "test";
    strcpy(file,fileprefix);
    const char* suffix;
    const char *mode;
    if (OutputFormat == "fa"){
      suffix = ".fa";
      mode = "wu";
    }
    if (OutputFormat == "fa.gz"){
      suffix = ".fa.gz";
      mode = "wb";
    }
    if (OutputFormat == "fq"){
      suffix = ".fq";
      mode = "wu";
    }
    if (OutputFormat == "fq.gz"){
      suffix = ".fa.gz";
      mode = "wb";
    }
    //strcat(file,suffix);
    // fprintf(stderr,"%s",file);
    int mt_cores = 1;
    samFile *SAMout = NULL;
    sam_hdr_t *SAMHeader;
    bam1_t *bam_file_chr = bam_init1();
    if (OutputFormat == "bam"){
      
      htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
      const char* filename2 = "test.bam";
      if ((SAMout = sam_open_format(filename2, "wb", fmt_hts)) == 0) {
        fprintf(stderr,"Error opening file for writing\n");
        //insert sam_mt 
        exit(0);
      }
      SAMHeader = sam_hdr_init();
      Header_func(fmt_hts,filename2,SAMout,SAMHeader,seq_ref,chr_total);
      fprintf(stderr,"\ncreating file works\n");
      free(fmt_hts);
    }

    BGZF *bgzf;
    if (OutputFormat == "fq"){
      const char* filename = "test.fq";
      bgzf = bgzf_open(filename,mode);
      bgzf_mt(bgzf,mt_cores,256);
      fprintf(stderr,"\t-> Number of cores for bgzf_mt: %d\n",mt_cores);  
    }
    
    // READ QUAL ARRAY
    double* Qual_freq_array = new double[6000];
    Qual_freq_array = Qual_array(Qual_freq_array,"/home/wql443/WP1/SimulAncient/Qual_profiles/Acc_freq1.txt");
  
    // FRAGMENT LENGTH CREATING ARRAY
    int* Frag_len = new int[4096];
    double* Frag_freq = new double[4096];
    int number;

    FragArray(number,Frag_len,Frag_freq,"Size_dist/Size_dist_sampling.txt");
    
    int maxsize = 5;
    //initialzie values that should be used for each thread
    for (int i = 0; i < nthreads; i++){
      //fprintf(stderr,"threads loop \n");
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;
      struct_for_threads[i].threadno = i;
      struct_for_threads[i].genome = genome_data;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;

      struct_for_threads[i].FragLen =Frag_len;
      struct_for_threads[i].FragFreq = Frag_freq;
      struct_for_threads[i].No_Len_Val = number; 
      //struct_for_threads[i].sizearray = sizearray;
      //struct_for_threads[i].SizeDist = SizeDist;

      struct_for_threads[i].Qualfreq = Qual_freq_array;
      struct_for_threads[i].reads = reads;

      struct_for_threads[i].OutputFormat = OutputFormat;
      struct_for_threads[i].bgzf = bgzf;
      struct_for_threads[i].SAMout = SAMout;
      struct_for_threads[i].SAMHeader = SAMHeader;
      //struct_for_threads[i].BAM = bam_file_chr;

      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      
      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (int*)malloc(sizeof(int) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      /*for(int jj=0;jj<chr_total+1;jj++)
      struct_for_threads[i].size_cumm[jj]= chr_size_cumm[jj];*/
      
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));

      struct_for_threads[i].l = 0;
      struct_for_threads[i].m = 20;
      for(int j=0; j<20;j++){
        struct_for_threads[i].list_of_reads[j]=bam_file_chr;
        //fprintf(stderr,"for loop inside threads\n");
      }
    }
    // free(catString);
    // delete[] sizearray;
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    // pthread_create( thread ID, struct for thread config, pointer to the function which new thread starts with, 
    // it returns a pointer to void data which is the 4th argument 
    fprintf(stderr,"Creating a bunch of threads\n"); 
    for (int i = 0; i < nthreads; i++){
      pthread_create(&mythreads[i],&attr,Fafq_thread_se_run,&struct_for_threads[i]);
    }
    fprintf(stderr,"Done Creating a bunch of threads\n");
    fflush(stderr);

    for (int i = 0; i < nthreads; i++)
    {  
      fprintf(stderr,"joing threads\n");fflush(stderr);
      // PRINT 
      pthread_join(mythreads[i],NULL);
      //fprintf(stderr, "\t[ANDET STED] walltime used for join =  %.2f sec\n", (float)(time(NULL) - t3));  
    }

    fprintf(stderr,"Header close\n");
    sam_hdr_destroy(SAMHeader);
    sam_close(SAMout);
    fprintf(stderr,"Sam close");
    
    if (OutputFormat == "fq"){bgzf_close(bgzf);}
    fprintf(stderr,"AFTER loop\n");
    //for(int i=0;i<nthreads;i++){delete struct_for_threads[i].fqresult_r1 -> s;} //ERROR SUMMARY: 9 errors from 5 contexts (suppressed: 0 from 0 )
    for(int i=0;i<nthreads;i++){
      fprintf(stderr,"inside loop 1 \n");
      free(struct_for_threads[i].fqresult_r1 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      //ks_release(struct_for_threads[i].fqresult_r1);
      fprintf(stderr,"inside loop 2 \n");
      delete struct_for_threads[i].fqresult_r1;
      fprintf(stderr,"inside loop 3 \n");
      //delete[] struct_for_threads[i].Qualfreq;
      //delete[] struct_for_threads[i].sizearray;
      //free(struct_for_threads[i].fqresult_r1); // ERROR SUMMARY: 8 errors from 4 contexts (suppressed: 0 from 0)
    }
    //delete[] sizearray;
    fprintf(stderr,"outside of final loop\n");
    delete[] Frag_len;
    delete[] Frag_freq;
    delete[] Qual_freq_array;
    fprintf(stderr,"delete qual_freq_array \n");
    free(genome_data);
    fprintf(stderr,"free genome data \n");
  }
  fprintf(stderr,"before return statement \n");
  return NULL;
}

// ------------------------------ //
int main(int argc,char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);
  // printTime(stderr);
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa chr1_3 chr1 hg19canon.fa  chr1_12   chr22 chr1_15 chr10_15  chr18_19  
  //const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  fprintf(stderr,"\t-> fasta load \n");
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);
  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  int Glob_seed = 1; //(int) time(NULL);
  int threads = 1;
  size_t No_reads = 10;
  fprintf(stderr,"\t-> Seed used: %d with %d threads\n",Glob_seed,threads);
  fprintf(stderr,"\t-> Number of simulated reads: %zd\n",No_reads);
  //char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  const char* Adapt_flag = "false";
  const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
  
  const char* OutputFormat = "bam";
  //fprintf(stderr,"Creating a bunch of threads\n");
  int Thread_specific_Read = static_cast<int>(No_reads/threads);

  Create_se_threads(seq_ref,threads,Glob_seed,Thread_specific_Read,Adapt_flag,Adapter_1,OutputFormat);
  
  fprintf(stderr,"Done creating a bunch of threads\n");
  //fflush(stderr);
  // free the calloc memory from fai_read
  //free(seq_ref);
  fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  return 0;  
}

// g++ NGSNGS_func.cpp atomic_fq.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
// valgrind --tool=memcheck --leak-check=full ./a.out