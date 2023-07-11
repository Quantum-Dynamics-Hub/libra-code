/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_state_tracking.cpp
  \brief The file implements the algorithms for state tracking
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


int step1(MATRIX& X){
    /**
    For each row of the matrix, find the smallest element and subtract it from every 
    element in its row. Go to Step 2. 
    */

    int i,j;
    int n = X.n_cols;

    for(i=0; i<n;i++){
        //complex<double> res;
        double val;
        int min_elt_indx;
        X.min_row_elt(i, val, min_elt_indx); 

        for(j=0;j<n;j++){
            X.add(i, j, -val);  
        }
    }

    return 2;
}

int step2(MATRIX& X, vector< vector<int> >& M, vector<int>& Rcov, vector<int>& Ccov){
    /**
    Find a zero (Z) in the resulting matrix. If there is no starred zero in its row or column, 
    star Z. Repeat for each element in the matrix. Go to Step 3. 
    */
    int i,j;
    int n = X.n_cols;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(fabs(X.get(i,j))<1e-10 && Rcov[i]==0 && Ccov[j]==0){
                M[i][j] = 1;
                Rcov[i] = 1;
                Ccov[j] = 1;
            }
        }
    }// for i

    /**
    Before we go on to Step 3, we uncover all rows and columns so that we can use the
    cover vectors to help us count the number of starred zeros
    */

    for(i=0; i<n; i++){
        Rcov[i] = 0;
        Ccov[i] = 0;
    }
    
    return 3;    
}


int step3(vector< vector<int> >& M, vector<int>& Ccov){
    /**
    Cover each column containing a starred zero.  
    If K columns are covered, the starred zeros describe a complete set of unique assignments.
    In this case, Go to DONE, otherwise, Go to Step 4.
    */

    int i,j;
    int n = M.size();

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(M[i][j]==1){  Ccov[j] = 1; }
        }// for j
    }// for i

    int colcount = 0;
    for(i=0;i<n;i++){
        colcount += Ccov[i];
    }

    int step = 4;
    if(colcount>=n){   step = 7;   }

    return step;

}


void find_a_zero(MATRIX& X, vector<int>& Rcov, vector<int>& Ccov, int& row, int& col){

    int n = X.n_cols;
    int done = 0; // 0 - false, 1 - true
    row = -1;
    col = -1;

    int r = 0;
    while(!done){

        int c = 0;
        while(!(c>=n || done)){
            if(fabs(X.get(r,c))<1e-20 && Rcov[r]==0 && Ccov[c]==0){
                row = r;
                col = c;
                done = 1;
            }

            c += 1;
        }
 
        r += 1;
        if(r>=n){
            done = 1; // True
        }
    }

}


int find_star_in_row(vector< vector<int> >& M, int row){

    int n = M.size();

    int col = -1;    
    for(int i=0;i<n;i++){
        if(M[row][i] == 1){   col = i;  }
    }

    return col;
}


int find_star_in_col(vector< vector<int> >& M, int col){

    int n = M.size();

    int row = -1;
    for(int i=0;i<n;i++){
        if(M[i][col] == 1){  row = i; }
    }

    return row;

}


int find_prime_in_row(vector< vector<int> >& M, int row){
    int n = M.size();

    int col = -1;
    for(int i=0;i<n;i++){
        if(M[row][i] == 2){  col = i; }
    }

    return col;

}


int step4(MATRIX& X, vector< vector<int> >& M, vector<int>& Rcov, vector<int>& Ccov, int& path_row_0, int& path_col_0){
    /**
    Find a noncovered zero and prime it. 
    If there is no starred zero in the row containing this primed zero, Go to Step 5. 
    Otherwise, cover this row and uncover the column containing the starred zero. 
    Continue in this manner until there are no uncovered zeros left. 
    Save the smallest uncovered value and Go to Step 6.
    */

    int done = 0; // 0 - False, 1 - True
    int row = -1;
    int col = -1;
    int step = 6;

    while(!done){
        find_a_zero(X,Rcov, Ccov, row, col);
        if(row==-1){
            done = 1;
            step = 6;
        }
        else{
            M[row][col] = 2;

            int col_ = find_star_in_row(M, row);

            if(col_>-1){
                col = col_;
                Rcov[row] = 1;
                Ccov[col] = 0;
            }
            else{
                done = 1;
                path_row_0 = row;
                path_col_0 = col;
                step = 5;
            }
        }
    }

    return step;

}


void augment_path(vector< vector<int> >& M, vector< vector<int> >& path){

    int sz = path.size();

    for(int i=0; i<sz;i++){ 
        vector<int> p = path[i];

        if(M[p[0]][p[1]]==1){
            M[p[0]][p[1]] = 0;
        }
        else{
            M[p[0]][p[1]]=1;
        }
    }

}


void clear_covers(vector<int>& Rcov, vector<int>& Ccov){

    int n = Rcov.size();

    for(int i=0;i<n;i++){
        Rcov[i] = 0;
        Ccov[i] = 0;
    }
}


void erase_primes(vector< vector<int> >& M){

    int n = M.size();

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(M[i][j]==2){
                M[i][j] = 0;
            }
        }
    }

}



int step5(vector< vector<int> >& M, int path_row_0, int path_col_0, vector<int>& Rcov, vector<int>& Ccov, vector< vector<int> >& path){
    /**
    Construct a series of alternating primed and starred zeros as follows.
    Let Z0 represent the uncovered primed zero found in Step 4.
    Let Z1 denote the starred zero in the column of Z0 (if any). 
    Let Z2 denote the primed zero in the row of Z1 (there will always be one).  
    Continue until the series terminates at a primed zero that has no starred zero in its column.  
    Unstar each starred zero of the series, star each primed zero of the series, erase all primes
    and uncover every line in the matrix. Return to Step 3. 
    */
 
    int done = 0; // 0 - False, 1 - True
    int r = -1;
    int c = -1;

    
    path.clear();

    vector<int> pi = vector<int>(2, 0);
    pi[0] = path_row_0; pi[1] = path_col_0;
    path.push_back(pi);

    int path_count = 1;

    while(!done){
     
        r = find_star_in_col(M, path[path_count-1][1]);
        if(r>-1){
            path_count += 1;
            pi[0] = r; pi[1] = path[path_count-2][1];
            path.push_back(pi);
        }
        else{
            done = 1;
        }
      
        if(!done){
            c = find_prime_in_row(M, path[path_count-1][0]);
            path_count += 1;
            pi[0] = path[path_count-2][0]; pi[1] = c;
            path.push_back(pi);
        }
    }

    augment_path(M, path);
    clear_covers(Rcov, Ccov);
    erase_primes(M);

    return 3;

}



double find_smallest(MATRIX& X, vector<int>& Rcov, vector<int>& Ccov){

    int i,j;
    int n = X.n_cols;
    double minval = 1.0e+25;
  
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(Rcov[i]==0 && Ccov[j]==0){
                if(minval>X.get(i,j)){
                    minval = X.get(i,j);
                }
            }
        }
    }

    return minval;

}



int step6(MATRIX& X, vector<int>& Rcov, vector<int>& Ccov){
    /**
    Add the value found in Step 4 to every element of each covered row, and subtract it from
    every element of each uncovered column.  Return to Step 4 without altering any stars, 
    primes, or covered lines. 
    */

    double minval = find_smallest(X, Rcov, Ccov);

    int i,j;
    int n = X.n_cols;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
    
            if(Rcov[i]==1){
                X.add(i,j, minval);
            }
            if(Ccov[j]==0){
                X.add(i,j, -minval);
            }
        }
    }

    return 4;

}


vector<vector<int> > compute_assignment(vector<vector<int> >& M){

    vector< vector<int> > res;
    int i,j;
    int n = M.size();

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){

            if(M[i][j] == 1){
                vector<int> p(2, 0);
                p[0] = i; p[1] = j;
                res.push_back(p);
            }
        }
    }

    return res;

}


void show_M(vector<vector<int> >& M){

    int i,j;
    int n = M.size();

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            cout<<M[i][j]<<"  ";
        }
        cout<<"\n";
    } 
}


vector<int> get_permutation(vector<vector<int> >& inp){
    /**
    Convert the list of lists into the permutation object
 
    inp - N x N matrix (vector of vectors)
    */ 

    int sz = inp.size();
    vector<int> perm_t(sz, 0);

    for(int i=0;i<sz;i++){
        vector<int> ra = inp[i];
        perm_t[ra[0]] = ra[1];  // for < alpha | alpha > this becomes a new value: perm_t = P_{n+1}
    }

    return perm_t;
}




vector<int> Munkres_Kuhn_minimize(MATRIX& _X, int verbosity){

    MATRIX X(_X);

    int i,j;
    int n = X.n_cols;

    vector<vector<int> > M(n, vector<int>(n, 0));
    vector<int> Rcov(n, 0);
    vector<int> Ccov(n, 0);

    int done = 0; // 0 - False, 1 - True
    int step = 1;
    int path_row_0, path_col_0; // = None, None

    vector<vector<int> > path; 
    vector<vector<int> > res; 

    int iter_ = 2;

    while(!done){

        if(verbosity > 0){
            cout<<"** Iteration "<<iter_<<" **\n"; 
        }
        iter_ += 1;


        if(step==1){
            if(verbosity > 0){  cout<<"Step 1\n"; }
            step = step1(X);
        }

        else if(step==2){
            if(verbosity > 0){  cout<<"Step 2\n"; }
            step = step2(X, M, Rcov, Ccov);
        }

        else if(step==3){
            if(verbosity > 0){  cout<<"Step 3\n"; }
            step = step3(M, Ccov);
        }

        else if(step==4){
            if(verbosity > 0){  cout<<"Step 4\n"; }
            step = step4(X, M, Rcov, Ccov, path_row_0, path_col_0);
        }

        else if(step==5){
            if(verbosity > 0){  cout<<"Step 5\n"; }
            step = step5(M, path_row_0, path_col_0, Rcov, Ccov, path);
        }

        else if(step==6){
            if(verbosity > 0){  cout<<"Step 6\n"; }
            step = step6(X, Rcov, Ccov);
        }

        else if(step==7){
            res = compute_assignment(M);
            done = 1;

            if(verbosity > 0){  
                cout<<"Done\n";                    
                for(i=0;i<res.size();i++){
                    cout<<res[i][0]<<"  "<<res[i][1]<<" \n";
                }
            }

        }

        if(verbosity > 0){
            cout<<"X = \n"; X.show_matrix();
            cout<<"M = \n"; show_M(M);
            cout<<"Ccov = \n";
            for(i=0;i<Ccov.size();i++){ cout<<Ccov[i]<<" "; } cout<<"\n";
            for(i=0;i<Rcov.size();i++){ cout<<Rcov[i]<<" "; } cout<<"\n";
            cout<<"============================\n";
        }

    }// while

    return get_permutation(res);

}



vector<int> Munkres_Kuhn_maximize(MATRIX& _X, int verbosity){
    /**
    Minimize the negative of the original matrix
 
    We also need to shift the negative matrix up rigidly, so all its 
    elements are > 0.0
    */

    MATRIX X(_X);

    int i,j;
    int n = X.n_cols;
    double maxval = 0.0;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(X.get(i,j)>maxval){
                maxval = X.get(i,j);
            }
        }
    }

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            X.scale(i,j, -1.0);
            X.add(i,j, maxval + 1e-5);
        }
    }

    return Munkres_Kuhn_minimize(X, verbosity);

}


}// namespace libdyn
}// liblibra


