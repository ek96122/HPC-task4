//
//  main.cpp
//  HPC q4
//
//  Created by Ernest on 24/3/2016.
//  Copyright © 2016 Kan Tsz Hang. All rights reserved.
//

#include <iostream>
#include <vector>
#include "TriMatrix.h"
#include <Accelerate/Accelerate.h>

using namespace std;

int main(){
    double L, T, alpha, dt, theta;
    int Nx;
    
    
    // get user input
    {using namespace std;
        cout << "L: ";
        cin >> L;
        
        cout << "T: ";
        cin >> T;
        
        cout << "alpha: ";
        cin >> alpha;
        
        cout << "Nx: ";
        cin >> Nx;
        
        cout << "dt: ";
        cin >> dt;
        
        cout<< "theta (0 for forward Euler)/(0.5 for Crank-Nicolson)/(1 for backward Euler) :";
        cin >> theta;
    }
    
    double dx = L/Nx;
    double nu = alpha*dt/dx/dx;
    double gamma0 = 0;
    double gamma1 = 0;
    
    
    //Defining all BASIC vectors used in the programme
    
    double x[Nx+1];
    double u0[Nx+1];
    double *u1;
    u1 = new double[Nx+1];
    double *u2;
    u2 = new double [Nx+1];
    
    //Defining x-vector
    for (int i=0;i<Nx+1;i++) {
        x[i]=dx*i; //equallly spaced elements in x-vector
    }
    
    //Defining a vector, which uses the initial condition u0, to hold the solution
    for (int i=0; i<Nx+1; i++) {
        u0[i]=x[i]/(1-x[i]);
    }
    
    //Transfer entries from u0 vector to the u1 vector created earlier
    
    u1[0]=gamma0;
    u1[Nx]=gamma1;
    
    for (int i=1; i<Nx; i++) {
        
        u1[i] = u0[i];
    }
    
    
    
    TriMatrix Left(Nx, theta*nu*-1);
    TriMatrix Right(Nx, (1-theta)*nu);
    Right.convertMatrix();
    
    
    //forward Euler time integration via matrix multiplication
    
    for (double i=0; i<(T-dt); i+=dt){
        
        u2=Left.inverse(Right.blasmultiply(u1));
        u1=u2;
        cout<<i<<endl;
    }
    
    cout<<endl;
    cout <<"time-integration output" <<endl;
    for (int i=0; i< Nx + 1; i++){
        cout<<u1[i] << endl;
    }

    
    
    

    
    
    return 0;
    
    
}
