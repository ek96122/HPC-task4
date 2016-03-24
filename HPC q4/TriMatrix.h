//
//  TriMatrix.h
//  HPC q4
//
//  Created by Ernest on 24/3/2016.
//  Copyright © 2016 Kan Tsz Hang. All rights reserved.
//

#ifndef TriMatrix_h
#define TriMatrix_h
#include <vector>
#include <iostream>
#include <Accelerate/Accelerate.h>


using namespace std;

class TriMatrix{
    
private:
    double *diag;
    double *u_diag;
    double *l_diag;
    double *A;
    int length;
    
public:
  
    TriMatrix(int Nx, double nu){
        l_diag = new double[Nx];
        diag = new double[Nx+1];
        u_diag = new double[Nx];
        
        for (int i = 0; i < Nx + 1; i++){
            diag[i] = 1-2*nu;
        }
        
        for (int i = 0; i < Nx; i++){
            l_diag[i] = nu;
            u_diag[i] = nu;
        }
        diag[0] = 1;
        diag[Nx] = 1;
        l_diag[Nx-1] = 0;
        u_diag[0] = 0;
        length = Nx + 1;
        
    }
    
    void convertMatrix(){
        A = new double [length*length];
        
        for (int i=0; i<length; i++) {
            A[i]=0;
        }
        
        for (int i=0; i<length; i++){
            A[i+i*length]=diag[i];
        }
    
        for (int i=0; i< length-1; i++){
            A[i+i*length+1]=l_diag[i];
        }
    
        for (int i=0; i<length-1; i++){
            A[i+(i+1)*length]=u_diag[i];
        }
    }
    
        
    
    
    //Create a function which can perform matrix multiplication for the tridiagonal matrix with any Nx1 matrix, where N=number of rows of the tridiagonal matrix
        double *blasmultiply(double *x){
        int alpha = 1;
        int beta = 0;
        double *r;
        
        r = new double[length];
    
        cblas_dgemv(CblasColMajor, CblasNoTrans, length, length, alpha, A, length, x, 1, beta, r, 1);
        
      
        return r;
    }
    
    
    double *inverse(double *x){
        int nrhs = 1;
        int info;
        double *dl, *d, *du;
        dl = new double[length-1];
        d = new double[length];
        du = new double[length-1];
        copy(l_diag, l_diag + length-1, dl);
        copy(diag, diag + length, d);
        copy(u_diag, u_diag + length-1, du);
        
        dgtsv_(&length, &nrhs, dl, d, du, x, &length, &info);
        
        return x;
        
    }
    
    
        
        //cout << endl << "U2: " << endl;
        //for(int i=0; i<D; i++) cout << (U2)[i] << ", ";
        //cout<<endl;
    
        
        

    
};



#endif /* TriMatrix_h */
