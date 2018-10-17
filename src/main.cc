#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>

#include "scDiffNet.h"
#include "mmio.h"


void print_usage (FILE * stream, int exit_code)
{
	fprintf (stream, "Usage: DECODE\n ");
    fprintf (stream,
		 "\t-h	--help\t\t\tDisplay this usage information. \n"
		 "\t-i	--input FILENAME\tfullpath of the input expression file in mm format (MANDATORY)\n"
		 "\t-t	--type {mm|table|csv}\tformat of input file (MANDATORY)\n"		 
		 "\t-g	--group COLFILE\t\tFile with each row describing columns of interest in the expression matrix\n"		 
		 "\t-n	--null COLFILE\t\tFile with each row describing null columns in the expression matrix (default=all columns minus \"group\")\n"		 
		 "\t-r	--rows ROWFILE\t\tRows of matrix to perform differential analysis (default = all rows)\n"		 
		 "\t-p	--thread_no INT\t\tnumber of threads (default = -1)\n"
		 "\t-w	--weighted DOUBLE\t\tZ-score threshold for filtering ubiquitous features (default = 3)\n"
		 "\t-o	--output PATH\t\tfullpath of the folder to store output files (default = ./results)\n"
	);	 
		
    exit (exit_code);
}


int main(int argc, char ** argv) {
  int next_option;
  const char *const short_options = "hi:t:p:o:zw:";
  const struct option long_options[] = {
		{"help",     0, NULL, 'h'},
		{"input",  	 1, NULL, 'i'},
		{"type",  	 1, NULL, 't'},
		{"thread_no", 	 1, NULL, 'p'},
		{"weighted", 	 1, NULL, 'w'},
		{"output", 	 1, NULL, 'o'},
		{NULL,       0, NULL,  0 }		
	};


	char output_path[1024] = "./results", input_file[1024] = "";
	bool weighted = true;
		
	int type = -1; // -1: unknown, 0: mm, 1: table, 2: csv
	int thread_no = -1;
    do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option) {
			case 'h':
	    		print_usage (stdout, 0);

			case 'i':
				strcpy(input_file, optarg);
	    		break;

			case 't':
				if(!strcmp(optarg, "mm"))
					type = 0;
				else if(!strcmp(optarg, "table"))
					type = 1;				
				else if(!strcmp(optarg, "csv"))
					type = 2;
				else {
					fprintf(stderr, "unknown file format\n");
					print_usage(stderr, -1);					
					return -1;
				}
	    		break;
   		
			case 'p':
				thread_no = atoi(optarg);
	    		break;			

			case 'w':
				weighted = atof(optarg) > 0;
	    		break;			

			case 'o':
				strcpy(output_path, optarg);
	    		break;	    		
		}
    } while (next_option != -1);


	sp_mat expression;
	if(!strcmp(input_file, "") || type == -1) {
		fprintf(stderr, "full path to the expression file and its type are mandatory arguments.\n");
		print_usage(stderr, -1);
		return -1;
	}
	if(group.n_elem == 0) {
		fprintf(stderr, "primary column group argument is mandatory (-g | --group).\n");
		print_usage(stderr, -1);
		return -1;
	}
	
	switch(type) {
		case 0:
			expression = read_from_mm(input_file);
			break;
			
		case 1:
			expression = read_from_table(input_file);
			break;
			
		case 2:
			expression = read_from_csv(input_file);
			break;
					
	}

	mat net = read_from_csv(input_file, weighted);

	field<mat> res = constructNet(expression, net, thread_no);
	
	return 0;
}
