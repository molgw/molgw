#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>

using namespace std;

struct basis{
 int l;
 double expon,coef;
};
vector<basis>basis_l_all;
vector<basis>basis_l;
vector<basis>basis_s;

int main(int argc, char *argv[])
{
 bool newprim;
 int ibasis,jbasis,iprim,jprim,nbasis_l,nprim,l,lp1,lm1,nprint;
 double coef,expon;
 vector<int>prims;
 vector<int>lls;
 if(argc==2)
 {
  string basis_file_in=argv[1];
  string basis_file_out=basis_file_in+"_rel";
  ifstream read_bas_file(basis_file_in);
  read_bas_file>>nbasis_l;
  for(ibasis=0;ibasis<nbasis_l;ibasis++)
  {
   read_bas_file>>nprim>>l;
   prims.push_back({nprim});
   lls.push_back({l});
   for(iprim=0;iprim<nprim;iprim++)
   {
    read_bas_file>>expon>>coef;
    basis_l_all.push_back({l,expon,coef});
    if(ibasis==0 && iprim==0)
    {
     basis_l.push_back({l,expon,coef});
    }
    else
    {
     newprim=true;
     for(jbasis=0;jbasis<basis_l.size();jbasis++)
     {
      if(basis_l[jbasis].l==l && basis_l[jbasis].expon==expon)
      {
       newprim=false;
       jbasis=basis_l.size();
      }
     }
     if(newprim)
     {
      basis_l.push_back({l,expon,coef});
     }
    }
   } 
  }
  read_bas_file.close();
  coef=1.0e0;
  for(ibasis=0;ibasis<basis_l.size();ibasis++)
  {
   if(basis_l[ibasis].l==0)
   {
    lp1=basis_l[ibasis].l+1;
    basis_s.push_back({lp1,basis_l[ibasis].expon,coef});
   }
   else
   {
    lp1=basis_l[ibasis].l+1;
    basis_s.push_back({lp1,basis_l[ibasis].expon,coef});
    lm1=basis_l[ibasis].l-1;
    basis_s.push_back({lm1,basis_l[ibasis].expon,coef});
   }
  }
  ofstream write_bas_file(basis_file_out);
  write_bas_file<<prims.size()+basis_s.size()<<endl;
  write_bas_file<<setprecision(10)<<fixed;
  ibasis=0;
  for(iprim=0;iprim<prims.size();iprim++)
  {
   write_bas_file<<prims[iprim]<<setw(6)<<lls[iprim]<<endl;
   for(jprim=0;jprim<prims[iprim];jprim++)
   {
    write_bas_file<<setw(20)<<basis_l_all[ibasis].expon<<setw(25)<<basis_l_all[ibasis].coef<<endl;
    ibasis++;
   } 
  }
  nprint=0;l=0;
  do
  {
   for(ibasis=0;ibasis<basis_s.size();ibasis++)
   {
    if(basis_s[ibasis].l==l)
    {
     write_bas_file<<1<<setw(6)<<basis_s[ibasis].l<<endl;
     write_bas_file<<setw(20)<<basis_s[ibasis].expon<<setw(25)<<basis_s[ibasis].coef<<endl;
     nprint++;
    }
   }
   l++;
  }while(nprint!=basis_s.size());
  write_bas_file.close();
 }
 return 0;
}
