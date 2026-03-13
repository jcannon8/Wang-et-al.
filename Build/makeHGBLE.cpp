// makeHGBLE.cpp : Output random HGBLE string based on sequence and 
// frequency table.
// John F. Cannon 5-15-13
#include "stdafx.h"
#include <malloc.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define BUFFER_SIZE 6000
#define MAX_PEPTIDE_LENGTH 200
typedef struct {	// HGBLEfrequency fields
	char Name[4];
	double pH,pG,pB,pL,pE;
} ResidueData;

int GetResIndex(const char* buffer, ResidueData* pRes);
char GetConfor(int i,ResidueData* pRes);
int ParseDatabase(char* filename, ResidueData* pRes);
int ParseSequence(char* filename,ResidueData* pRes,int* pResIndex);

int main(int argc, char* argv[]) {
	if(argc<3) {	// The first two parameters are manditory.
		fprintf(stdout,"%s Output random HGBLE string using sequence and table.\n",
			argv[0]);
		fprintf(stdout,"USAGE:\n%s <protein sequence> <frequency database>\n",
			argv[0]);
		fprintf(stdout,"protein sequence is in the form of a LEaP source.\n");
		return 0;
	}
	// Parse the database file.
	ResidueData* pRes=(ResidueData*)malloc(21*sizeof(ResidueData));
	ParseDatabase(argv[2],pRes);
	// Parse sequence
	int* pResIndex=(int*)malloc(MAX_PEPTIDE_LENGTH*sizeof(int));
	int numResidues=ParseSequence(argv[1],pRes,pResIndex);
	// Now create random conformer model.
	char string[MAX_PEPTIDE_LENGTH];
	srand( (unsigned)time(NULL));	// Seed random number generator.
	string[0]=0;	// Empty buffer
	int ProIndex=GetResIndex("PRO",pRes);
	int j;
	for(j=0;j<numResidues;j++) {
		if(pResIndex[j]==ProIndex) string[j]='E';	// Leave prolines extended or what ever LEaP uses.
		else string[j]=GetConfor(pResIndex[j],pRes);
	}
	string[j]=0;	// NULL terminate 
	puts(string);	// Output string
}

int GetResIndex(const char* buffer, ResidueData* pRes) {
	// Return the database index for the residue name in buffer.
	int j;
	for(j=0;j<21;j++) {
		if(!strncmp(buffer,pRes[j].Name,3)) return j;
	}
	// If no match, return 21.
	return j;
}
char GetConfor(int i,ResidueData* pRes) {
	// return a HGBLE character for a residue of index i.
	// 0 to RAND_MAX is the range for rand(), need to bring that to 0 to 1.
	double Q=rand()/(double) RAND_MAX;
	// pH,pG,pB,pL,pE;
	if(Q<=pRes[i].pH) return 'H';
	if(Q<=pRes[i].pH+pRes[i].pG) return 'G';
	if(Q<=pRes[i].pH+pRes[i].pG+pRes[i].pB) return 'B';
	if(Q<=pRes[i].pH+pRes[i].pG+pRes[i].pB+pRes[i].pL) return 'L';
	return 'E';
}
int ParseDatabase(char* filename, ResidueData* pRes) {
	// Fill the pRes structure with data from frequency database in file, filename.		
	FILE* fpDatabase=fopen(filename,"r");
	if(!fpDatabase) {
		fprintf(stderr,"Problem opening database.\n");
		exit(0);
	}
	char buffer[BUFFER_SIZE];
	fgets(buffer,BUFFER_SIZE,fpDatabase);	// First line has the header information
	int i=0;
	char sH[6],sG[6],sB[6],sL[6],sE[6];
	while(6==fscanf(fpDatabase,"%3s %s %s %s %s %s",
		pRes[i].Name,sH,sG,sB,sL,sE)){
		pRes[i].pH=atof(sH);
		pRes[i].pG=atof(sG);
		pRes[i].pB=atof(sB);
		pRes[i].pL=atof(sL);
		pRes[i].pE=atof(sE);
		i++;
		if(i>21) {
			fprintf(stderr,"Too much data! Expected only 21 residues, found %d.\n",i);
			exit(0);
		}
	}
	fclose(fpDatabase);		// Done with file.
	if(i==21) return 1;
	fprintf(stderr,"Missing data only %d residues, expected 21.\n",i);
	exit(0); return 0;
}
int ParseSequence(char* filename,ResidueData* pRes,int* pResIndex) {
	// Parse the peptide sequence into a string of indices in pResIndex.
	FILE* fpSeq=fopen(filename,"r");
	// The protein sequence can be in the form of a LEaP source file like
	// hI2seg = sequence { NTYR HIS PRO ALA ASP LYS ASP TYR GLY LEU CASP }
	// Note the "{" and "}" which delimit the ends.
	// The terminal residues also have N or C prefixes. 
	// They are ignored and will end up extended.'
	char buffer[BUFFER_SIZE];
	int j=0,numResidues=0;
	while('{' != fgetc(fpSeq));	// Advance file pointer to after open brace.
	fscanf(fpSeq,"%s",buffer);	// Read first residue and ignore
	// Encode the sequence as a list of database indices in pResIndex of length numResidues.
	while(fscanf(fpSeq,"%s",buffer)) {
		if(!strcmp(buffer,"}")) break;	// End upon closing brace
		// Determine the database index to use for this residue.
		pResIndex[j]=GetResIndex(buffer,pRes);
		if(pResIndex[j]==21) break;
		// If residue name not in database, assume last residue.
		// Don't use HIE, HID, HIP becauase they aren't in the database. 
		j++;
	}
	numResidues=j;
	fclose(fpSeq);	// Done with file.
	// Residues that precede proline are XPR in the database.
	// Reassign the residue index for those that precede proline.
	int ProIndex=GetResIndex("PRO",pRes);
	int XprIndex=GetResIndex("XPR",pRes);
	for(j=0;j<numResidues;j++) {
		if(pResIndex[j]==ProIndex && j>0) 
			// At a proline, the preceeding residue is a XPR.
			pResIndex[j-1]=XprIndex;
	}
	return numResidues;
}

