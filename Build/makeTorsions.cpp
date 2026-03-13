// makeTorsions.cpp: Write tleap impose statements based on HGBLE string.
// John F. Cannon 5-15-13
#include "stdafx.h"
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#define BUFFER_SIZE 6000

int writeScript(char* imposeScript,char* string);

int main(int argc, char* argv[]) {
	if(argc<2) {	// The first two parameters are manditory.
		fprintf(stdout,"%s writes tleap impose statements based on HGBLE strings.\n",
			argv[0]);
		fprintf(stdout,"USAGE:\n%s <HGBLE string>\n",argv[0]);
		return 0;
	}
	char imposeScript[10*BUFFER_SIZE];
	writeScript(imposeScript,argv[1]);
//	puts(imposeScript);		// Output imposeScript
	fprintf(stdout,"%s",imposeScript); // Output imposeScript (no newline)
}
int writeScript(char* imposeScript,char* string) {
	// Construct LEaP impose commands based on HGBLEstring.
	char buf1[BUFFER_SIZE],buf2[BUFFER_SIZE];
	buf2[0]=0;	// Empty buffer
	int len=strlen(string);
	char begining[]="impose hI2seg {";
	// These regions of the Ramachandran plot define local structures.
	// The impose commands alter psi and phi appropriately.
	char alpha[]="} {{ C N CA C -57.0 }{ N CA C N -47.0 } }\n";
	char Ghelix[]="} {{ C N CA C -49.0 }{ N CA C N -26.0 } }\n";
	char beta[]="} {{ C N CA C -119.0 }{ N CA C N 113.0 } }\n";
	char loop[]="} {{ C N CA C 57.0 }{ N CA C N 47.0 } }\n";
	// Let LEaP leave other residues with the default extended
	// First look for H's
	int i;
	for(i=0;i<len;i++) {
		if(string[i]=='H') {
			sprintf(buf1,"%d ",i+2);	// The HGBLE string starts at residue 2.
			strcat(buf2,buf1);
		}
	}
	if(strlen(buf2)) {
		strcat(imposeScript,begining);
		strcat(imposeScript,buf2);
		strcat(imposeScript,alpha);
	}
	buf2[0]=0;	// Empty buffer
	// Next look for G's
	for(i=0;i<len;i++) {
		if(string[i]=='G') {
			sprintf(buf1,"%d ",i+2);	// The HGBLE string starts at residue 2.
			strcat(buf2,buf1);
		}
	}
	if(strlen(buf2)) {
		strcat(imposeScript,begining);
		strcat(imposeScript,buf2);
		strcat(imposeScript,Ghelix);
	}	
	buf2[0]=0;	// Empty buffer
	// Next look for B's
	for(i=0;i<len;i++) {
		if(string[i]=='B') {
			sprintf(buf1,"%d ",i+2);	// The HGBLE string starts at residue 2.
			strcat(buf2,buf1);
		}
	}
	if(strlen(buf2)) {
		strcat(imposeScript,begining);
		strcat(imposeScript,buf2);
		strcat(imposeScript,beta);
	}	
	buf2[0]=0;	// Empty buffer
	// Finally, look for L's
	for(i=0;i<len;i++) {
		if(string[i]=='L') {
			sprintf(buf1,"%d ",i+2);	// The HGBLE string starts at residue 2.
			strcat(buf2,buf1);
		}
	}
	if(strlen(buf2)) {
		strcat(imposeScript,begining);
		strcat(imposeScript,buf2);
		strcat(imposeScript,loop);
	}	
	return 1;
}

