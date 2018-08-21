#include"nlib.h"

double sq(double i){
   return i*i;
}

int renumnishi(  Inp_nishi inp1 ){
   puts("Start renumnishi");
   
   string pdbname = inp1.read("INPUTPDB1");
   string outfile1 = inp1.read("OUTPUTFILE1");
   int start_res = atoi(inp1.read("START_RES").c_str());

   FILE *fout1;
   if((fout1 = fopen(outfile1.c_str(),"w")) == NULL ){
      printf("cannot open output file: %s\n",outfile1.c_str());
      exit(1);
   }

   pdb_nishi* pdb1;
   pdb1 = new pdb_nishi(pdbname.c_str());
   cout<<"TOTAL ATOM = "<< pdb1->total_atom<<endl;
   cout<<"TOTAL RESIDUE = "<< pdb1->total_residue<<endl;
   
   int re_resnum = start_res;
   for(unsigned int i=0;i<pdb1->total_atom;i++){
      if( i < pdb1->total_atom - 1 && pdb1->rnum[i] != pdb1->rnum[i+1] ){
         pdb1->rnum[i] = re_resnum;
         re_resnum++;
      }
      else{
         pdb1->rnum[i] = re_resnum;
      }
   }

   pdb1->write_pdb(outfile1.c_str());

   fclose(fout1);
   return 0;
}
