#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <libgen.h>
#include <assert.h>




#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int get_nlines(FILE* fp)
{
    int nlines = 0;

    int ch;
    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == '\n') {
            nlines++;
        }
    }

    fseek(fp, 0, SEEK_SET);
    return nlines;
}

static char* get_value(char* s)
{
    char* pch;
    pch = strtok(s, " ");
    if (pch != NULL) {
        pch = strtok(NULL, " ");
        if (pch != NULL) {
            return strdup(pch); // need to free() after use
        } else {
            return NULL;
        }
    } else {
        return NULL;
    }
}

int main()
{
    FILE* fp = fopen("../mat", "r");
    if (fp != NULL) {
        int nlines = get_nlines(fp);

        printf("nlines: %d\n", nlines);

        // make values..
        char** values = calloc(nlines + 1, sizeof(char*)); // +1 for a sentinel
        if (values != NULL) {
            char line[1024];
            int idx = 0;
            while (fgets(line, sizeof(line), fp) != NULL) {
                values[idx] = get_value(line);
                idx++;
            }
        }

        // use values..
        char** p = &values[0];
        while (*p != NULL) {
            printf("%s\n", *p);
            p++;
        }

        // clean values..
        p = &values[0];
        while (*p != NULL) {
            free(*p);
            p++;
        }

        fclose(fp);
    } else {
        perror("../mat");
    }
}

/*



double** loadMattrix(int sss)
{   
    printf("here");
    FILE *ifile;
    
	ifile = fopen ("../mat", "r");
	double **D, val;
	int i, j, zsize;
	char *line = NULL;
	char *tok;
	char **tokens;
	
	// Read first line to get matrix size
	zsize = 15;
		
	D = (double **) calloc (zsize, sizeof (double *));
	for (i=0; i<zsize; i++)
	{
		D[i] = (double *) calloc (zsize, sizeof (double));
		memset (D[i], 0, (unsigned long) zsize * sizeof(double));
	}
    
	// Read matrix
	for (i=0; i<zsize; i++)
	{
		// Read line, remove trailing EOL, replace tabulation by blank space
		// line length = MAX_NAME_LENGTH + 16 blank spaces + ( N taxa * ([DECIMAL_DIG] digits + 2 blank spaces) )
		line = getLine (ifile, (40 + 16 + (zsize * (DECIMAL_DIG+2))));
		// Remove trailing EOL & replace tabulation by blank space
		line = str_replace ( line, '\n', "" );
		line = str_replace ( line, '\r', "" );
		line = str_replace ( line, '\t', " " );
		if (NULL != line)
		//if (NULL != (line = str_replace ( str_replace ( str_replace ( getLine (ifile, line, (MAX_NAME_LENGTH + 16 + (zsize * (DECIMAL_DIG+2)))), '\n', " " ), '\r', " " ), '\t', " " )))
		{
			// Split line on blank spaces
			tokens = str_split (line, ' ');
			if (tokens)
			{
				// Read each token
				for (j=0; j<=zsize; j++)
				{

                    tok = tokens[j];
						// First token is the sequence name
                    val = atof (tok);
                    D[j-1][i] = val;

                    free (tokens[j]);
	
				}
			}
			free(tokens);
		}
	}
	
	if (NULL != line)
		free (line);
	
	return D;
}

void mandi(){
    printf("mandi");
}

int main (int argc, char **argv) {
    printf("ok\n");
    // Options *options;
    // options = chooseSettings (argc, argv);
    double** D;
    printf("ok\n");
    D = loadMattrix(3);
    // D = loadM (options->fpI_data_file, &numSpecies, species);
    printf("ok\n");
    // run(argc, argv);
}

*/

