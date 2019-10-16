//matrix-coupledlog 23 June 2015
//this program should make a matrix of the coupled logistic equation transistion rates
//then "invert it" to solve the backwards Kolmogorov equation and the Kramers solution
//double = bad, double = better; but maybe this is the limit? - overflow errors!
//K from 2 to 50, a from 0.1 to 1.0 took 2600s
//K from 2^0 to 2^7, buffer 4, a 0.0 to 1.0 took 3934s (mostly from the last set); problem with 2^8, bad allocation?
//K from 2^0 to 2^7, buffer 5, a 0.0 to 1.0 took 9539s (~90% from the last set)
//K = 2^8 simply doesn't work, crashes after 10 minutes (without a useful error message)
//update 2016.04.14
//kay = double(8.+8.*jay) up to jay=16, so up to K=136 but it crashed but anyway gave negatives (ie didn't work)

#include <iostream>
#include <fstream>
#include <string>
//#include <typeinfo> //for: typeid - unnecessary?
#include <cmath> //for: pow
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include <complex>//unnecessary?
//#include <algorithm>//for: sort

using namespace std;
//using namespace arma;
using namespace Eigen;


//define birth and death rates
double bn (int n,int m, double params[6]) {//should probably include "[6]" in all of them so no/less junk gets in
    return params[0]*double(n);} // ! should all the inputs be doubles? or at least the output specifically?
double bm (int n,int m, double params[]) {
    return params[1]*m;}
double dn (int n,int m, double params[]) {
    return params[0]*n*(n + params[4]*m)/params[2];}
double dm (int n,int m, double params[]) {
    return params[1]*m*(params[5]*n + m)/params[3];}
double gg (int n,int m, double params[]) {
    return bn(n,m,params) + bm(n,m,params) + dn(n,m,params) + dm(n,m,params);}
double gn (int n,int m, double params[]) {
    return bm(n,m,params) + dn(n,m,params) + dm(n,m,params);}
double gm (int n,int m, double params[]) {
    return bn(n,m,params) + dn(n,m,params) + dm(n,m,params);}
double gnm (int n,int m, double params[]) {
    return dn(n,m,params) + dm(n,m,params);}

int vectorarg (int n,int m, int sumsize) {
    return int(m*sumsize-sumsize+n-1);}
double initconds (double params[6]) {//returns one double of K/(1+a), that is, where the fixed point is
    double tempini = params[2]/(1.+params[4]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}
double initcondx (double params[6]) {//returns x f.p.
    double tempini = (params[2]-params[4]*params[3])/(1.-params[4]*params[5]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}
double initcondy (double params[6]) {//returns y f.p.
    double tempini = (params[3]-params[5]*params[2])/(1.-params[4]*params[5]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}
double initcondx2skewed (double params[6]) {//returns x f.p.
    double tempini = params[2]/(1.+2.*params[5]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}
double initcondy2skewed (double params[6]) {//returns y f.p.
    double tempini = params[3]/(1.+2.*params[5]);
    if (tempini>1) {return tempini;}
    else {return 1.0;}
}

//double fudgefactor = 1.0e13; //included to prevent overflow error - grow the matrix to shrink its inverse

SparseMatrix<double> makesparsetrix(int interspec, double params[6]){    //makes the sparse matrix
//    could also use Triplets (http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#TutorialSparseFilling)
    int msize = interspec*interspec;
    SparseMatrix<double> mat(msize,msize);
    cout << "[square] matrix size = " << mat.rows() << endl;

    mat.reserve(VectorXi::Constant(msize,5));//each column will have at most five elements
    for (int j=0; j<interspec-1; j++){
          for (int l=0; l<interspec-1; l++){
                    int x=interspec*j+l; int y=interspec*j+l;
                    mat.insert(x,y)=-gg(l+1,j+1,params);
                    mat.insert(x+1,y)=bn(l+1,j+1,params);
                    mat.insert(x,y+1)=dn(l+2,j+1,params);
                    mat.insert(x+interspec,y)=bm(l+1,j+1,params);
                    mat.insert(x,y+interspec)=dm(l+1,j+2,params);
    }};
    for (int j=0; j<interspec-1; j++){
                int x=interspec*j+interspec-1;
                int y=interspec*(interspec-1)+j;
                mat.insert(x,x)=-gn(interspec,j+1,params);
                mat.insert(x+interspec,x)=bm(interspec,j+1,params);
                mat.insert(x,x+interspec)=dm(interspec,j+2,params);
                mat.insert(y,y)=-gm(j+1,interspec,params);
    };
    mat.insert(msize-1,msize-1)=-gnm(interspec,interspec,params);

    mat.makeCompressed();

//    mat = mat*fudgefactor;

    return mat;
}

// I guess I need to do K2 from 8 to 64, K1 from 80 to 640, and a21 from .1 to .9, using a12=1
// note also to do the loop aratio from .1 to .9 outermost, Kratio from aratio+.1 to .9, with new files written in each of these
int main()
{
    int buffer = 3; //how many times bigger than K is the artificial reflecting boundary (3 is okay, 4 is better, 5 is enough+)
    double rconst = 1.0;
    cout << "let us do " << (172-8)/8 << " of these..." << endl;

    for(double aaa=0.0; aaa<1.01; aaa+=0.1){
//     for(double jay=aaa+0.1;jay<0.91;jay+=0.1){
     for(double jay=0.5;jay<0.51;jay+=0.1){
      //define a file to write to
//      char numstr[round(10*jay)]; // enough to hold all numbers up to 64-bits
//      string result = "kramers-nosym-20180222-K2overK1is"+std::to_string(round(10*jay))+"a12is10a21is"+std::to_string(round(10*aaa))+".txt"
//      result = name + itoa(round(10*jay), numstr, 10);
      std::stringstream ss;
      ss << "kramers-nosym-20190730-K2overK1is" << round(10*jay) << "a12overa21is4a21is" << round(10*aaa) << ".txt";
      //ss << "kramers-nosym-20181213-K2overK1is1a12is0a21is" << round(10*aaa) << ".txt";
      //ofstream filer; filer.open("kramers-nosym-20180222-K1is%iK2is%ia12is10a21is%i.txt",round(kay1),k2,round(10*a21));
//      ofstream filer; filer.open("kramers-nosym-20180222-K1is"+to_string(round(kay1))+"K2is"+to_string(k2)+"a12is20a21is"+to_string(round(10*a21))+".txt");
      //ofstream filer; filer.open("kramers-nosym-20180222-K2overK1is"+std::to_string(round(10*jay))+"a12is10a21is"+std::to_string(round(10*aaa))+".txt");
      //ofstream filer; filer.open("kramers-nosym-20180222-K1is"+str(int(kay1))+"K2is"+str(k2)+"a12is20a21is"+str(int(10*a21))+".txt");
      ofstream filer; filer.open(ss.str());
//      ofstream filer; filer.open("temp.txt");
      for(int k2=8;k2<172+1;k2+=8){//got bad allocation error for 64,60 -> ran out of memory; 48 worked in 2 hours
        //define constants
        double a21 = aaa;
        double a12 = 0.5;//4.*a21;
        double kay2 = double(k2);
        double kay1 = double(k2);//int(k2/jay));
        double constants[6] = {rconst, rconst, kay1, kay2, a12, a21};
        int sumsize = round(constants[2])*(buffer);
        int matsize = sumsize*sumsize;//or could use pow

        //make and fill matrix
        SparseMatrix<double> mat = makesparsetrix(sumsize,constants);
//        SparseMatrix<double> matT = mat.transpose();
        mat.makeCompressed();
//        matT.makeCompressed();

//        //solving Mx=b, really M_b T=-1, so define vectors - don't forget to transpose M!
//        VectorXd tims;
//        VectorXd negonesvec = VectorXd::Constant(matsize,1,double(-1.0));
//        VectorXd fluxout = VectorXd::Zero(matsize,1);
//        for(int eye=0;eye<matsize;eye++){fluxout(eye)+=dn(eye,1,constants);fluxout(eye*sumsize)+=dm(1,eye,constants);};
        //or possibly we need to sum all of one row of M^-l //but this may need transposes, depending on whether I want row or column //Kramer
        VectorXd sumthese;//Kramer
        VectorXd pickaline = VectorXd::Zero(matsize,1);//Kramer
//        int startat = int(initconds(constants)); pickaline(vectorarg(startat,startat,sumsize))=1.0;//Kramer
        pickaline(vectorarg(int(initcondx2skewed(constants)),int(initcondy2skewed(constants)),sumsize))=1.0;//Kramer
//        pickaline(vectorarg(kay,kay,sumsize))=1.0;//Kramer

        //try solving it
        SparseLU<SparseMatrix<double>/**/, COLAMDOrdering<int>/**/ > solver;//not clear to me what the ordering is or why it's important
        solver.analyzePattern(mat);//this is also not clear to me why it's used, and/or if it's exclusive with the somment above
        solver.factorize(mat);
        sumthese = solver.solve(pickaline);//Kramer
//        solver.analyzePattern(mat.transpose());//this is also not clear to me why it's used, and/or if it's exclusive with the somment above
//        solver.factorize(mat.transpose());
//        tims = solver.solve(negonesvec);

//        SparseLU<SparseMatrix<double>/**/, COLAMDOrdering<int>/**/ > solver2;
//        solver2.analyzePattern(matT);
//        solver2.factorize(matT);
//        tims = -solver2.solve(negonesvec);

        //let's see
//        cout << "For K const " << kay << " and a " << aaa << " you get an extinction time of " << tims(int(kay*sumsize-sumsize+kay-1)) << endl;
//        cout << "For K const " << kay << " and a " << aaa << " you get an extinction time of " << -sumthese.sum() << " times " << fudgefactor << endl;
//        cout << endl;
        filer << k2 << "\t" << -sumthese.sum() << "\n";
//        filer << aaa << "\t" << jay << "\t" << -sumthese.sum() << "\n";
      }
      filer.close();
     }
    }

//    getchar();//to press enter - unnecessary
    return 0;
}

