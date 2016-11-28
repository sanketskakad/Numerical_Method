/*********************************************************************************************
  Author : Sanket Sanjay Kakad
  
  Information : This Library class will support following methods to solve roots
                
                1. Bisection Method
                2. Regula Falsi Method
                3. Newton Rhapson Method
**********************************************************************************************/
//#include "Numerical_Method.h"
#include <iostream>
/* MAX Define Maximum Degree of equation*/
#define MAX 100
/* Max_Bia_iteration Define Iteration by Bisection Method*/
#define Max_Bia_iteration 10
/* Max_NR_iteration Define Iteration by Newton Rhapson Method*/
#define Max_NR_iteration 10
/* Max_Bia_RFmethod Define Iteration by Regula Falsi Method*/
#define Max_Bia_RFmethod 10

using namespace std;

class Equation{
  public:
  /* Constructor */
  Equation()
  {
    int i;
    Len = 0;
    for(i=0; i< MAX ;i++)
    {
      coefficient[i] = 0;
    }
  }
  
  /**********************************************************************
  * Method Use : This method used for collecting coefficient of Equation*
  ***********************************************************************/
  void Enter_Coeff(void)
  {
    int i;
    cout <<"Enter Highest degree of Equation :";
    cin >> Len;
    cout <<endl;
    if(Len < MAX)
    {
      for(i=0; i<= Len ;i++)
      {
        cout << "Enter Coefficient Of Order " << i <<" :";
        cin >>coefficient[i];
      }
    }
    else
    {
      cout << endl <<"Out Of Range "<< endl;
    }
  }
  
  /**********************************************************************
  * Method Use : This method used for Printing coefficient of Equation  *
  ***********************************************************************/
  void Print_Coef(void)
  {
    int i;
    for(i=0; i<= Len ;i++)
    {
      cout << "Coefficient of order"<<i <<" is " << coefficient[i] <<endl;
    }
  }
  
  /**********************************************************************
  * Method Use : This method used for solving equation by Newton Rhapson *
  *              Method                                                  *
  ***********************************************************************/
  
  void NR_method()
  {
    double E , En;
    double EAns = 0, EAnsD = 0;
    int i, j;
    cout <<"Enter Fist Estimate :";
    cin >> E;
    for(j=0;j < Max_NR_iteration; j++)
    {
      EAns= 0;
      EAnsD = 0;
      for(i=0; i <= Len ;i++)
      {
        EAns += coefficient[i] * pow1(E,i)  ;
        EAnsD += i*coefficient[i] * pow1(E,(i-1));
      }
      if(EAnsD != 0 )
      {
        En = E - (EAns / EAnsD);
        cout << "At iteration "<<j << "\tF(X) = "<< EAns << "\tF'(X) = "<< EAnsD << "\tX"<<j<<" = "<< En << endl;
        E= En;
      }
      else
      {
        cout << "NR cant process further";
        break;
      }
    }
  }
  
  /**********************************************************************
  * Method Use : This method used for solving equation by Regula Falsi  *
  *              Method                                                 *
  ***********************************************************************/
  
  void RF_method ()
  {
    double L , U ;
    
    cout <<"Enter Lower Range :";
    cin >> L;
    
    cout <<"Enter Upper Range :";
    cin >> U;
    
    I_call_RFmethod(U,L,0);

  }
  
  /**********************************************************************
  * Method Use : This method used for solving equation by Bisection     *
  *              Method                                                 *
  ***********************************************************************/
  
  void BisectionMethod ()
  {
    double L , U ;
    
    cout <<"Enter Lower Range :";
    cin >> L;
    
    cout <<"Enter Upper Range :";
    cin >> U;

    I_call_Bisection(U,L, 0);
    
  }
  
  private:
  
  /**************************************************************************
  * Method Use : This Private method used for solving power related problem  *
  ***************************************************************************/
  
  double pow1(double a, int b)
  {
    double temp=1; 

    for(b;b>0;b--)
    {
        temp *= a;
    }
    return temp;
  }
  
  /**************************************************************************
  * Method Use : This Private method used for solving Regula Falsi method  *
  * Call Type : Recursive                                                  *
  ***************************************************************************/
  
  void I_call_RFmethod(double U , double L,int j)
  {
    int i;
    double UAns = 0;
    double LAns = 0;
    double MAns = 0;
    double M;
    if( j < Max_Bia_RFmethod)
    {
      for(i=0; i <= Len ;i++)
      { 
        UAns += coefficient[i] * pow1(U,i) ; 
      }
      
      for(i=0; i <= Len ;i++)
      {
        LAns += coefficient[i] * pow1(L,i)  ;
      }
      M = U - (((U-L)/(UAns-LAns))*UAns);
      
      for(i=0; i <= Len ;i++)
      {
        MAns += coefficient[i] * pow1(M,i)  ;
      }
      cout<<"A =" <<U <<"\tB =" << L << "\tM ="<< M << "\tf(A) =" << UAns << "\tf(B) =" << LAns<< "\tf(M) ="<< MAns <<endl;
      if(UAns*LAns < 0)
      {
        cout << "Iteration " << j << "Value is " << M;
        if(UAns*MAns < 0.000)
        {
          I_call_RFmethod(U,M,++j);
        }
        else if(LAns*MAns < 0.000)
        {
          I_call_RFmethod(M,L,++j);
        }
        else
        {
          cout << "Root Is "<<M;
        }
      }
      else
      {
        cout << endl << "For given interval either root is absent or there are even roots are present" <<endl;
      }
    }
  }
    

  /**************************************************************************
  * Method Use : This Private method used for solving Bisection method      *
  * Call Type : Recursive                                                   *
  ***************************************************************************/
  void I_call_Bisection(double U , double L,int j)
  {
      int i;
      double UAns = 0;
      double LAns = 0;
      double MAns = 0;
      double M;
    
    if( j <Max_Bia_iteration)
    {
      
      for(i=0; i <= Len ;i++)
      { 
        UAns += coefficient[i] * pow1(U,i) ; 
      }
      
      for(i=0; i <= Len ;i++)
      {
        LAns += coefficient[i] * pow1(L,i)  ;
      }

      M = 0.5 * ( U + L );
      
      for(i=0; i <= Len ;i++)
      {
        MAns += coefficient[i] * pow1(M,i)  ;
      }
      cout<<"\nA =" <<U <<"\tB =" << L << "\tM ="<< M << "\tf(A) =" << UAns << "\tf(B) =" << LAns<< "\tf(M) ="<< MAns <<endl;
      if(UAns*LAns < 0)
      {
        cout << "\nIteration " << j << " Value is " << M <<endl;
        if(UAns*MAns < 0.000)
        {
          I_call_Bisection(U,M,++j);
        }
        else if(LAns*MAns < 0)
        {
          I_call_Bisection(M,L,++j);
        }
        else
        {
          cout << "Root Is "<<M;
        }
      }
      else
      {
        cout << endl << "For given interval either root is absent or there are even roots are present" <<endl;
      }
    }
    else
    {
        
    }
  }
  /******************************************************************************************************
  **************************            PRIVATE VARIABLE             ************************************
  ******************************************************************************************************/
  
  int Len;
  double coefficient[MAX];
};

int main(void)
{
  Equation a;
  a.Enter_Coeff();
  cout << "\n**************************************************************\n";
  a.Print_Coef();
  cout << "\n**************************************************************\n\n";
  cout << "\n\n\n By using Bisection Method :\n";
  a.BisectionMethod();
  cout << "\n**************************************************************\n\n";
  cout << "\n\n\n By using Newton Rhapson Method :\n";
  a.NR_method();
  cout << "\n**************************************************************\n\n";
  cout << "\n\n\n By using Regula Faulsi Method : \n";
  a.RF_method();
  
  return 0;
}




