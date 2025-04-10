/*********************************************************************************
* Copyright (C) 2025 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file sHamiltonian.cpp
  \brief The file implements the key methods of the sHamiltonian class
    
*/

#if defined(USING_PCH)
#include "../pch.h"
#else
#include <stdlib.h>
#endif 

//#include <torch/torch.h>
#include "sHamiltonian.h"
#include <torch/csrc/api/include/torch/python.h>

#include "../math_meigen/mEigen.h"

namespace bp = boost::python;

/// liblibra namespace
namespace liblibra{

using namespace libmeigen;

/// libnhamiltonian namespace 
namespace libshamiltonian{



bp::object tensor_to_python(const torch::Tensor& tensor) {
    // Convert torch::Tensor to a NumPy array, then to a Python object
    PyObject* array = THPVariable_Wrap(tensor);
    return bp::object(bp::handle<>(array));
}

torch::Tensor python_to_tensor(const bp::object& obj) {
    PyObject* py_obj = obj.ptr();
    if (!THPVariable_Check(py_obj)) {
        throw std::runtime_error("Expected a PyTorch tensor.");
    }
    return THPVariable_Unpack(py_obj);
}




sHamiltonian::sHamiltonian(const sHamiltonian& src){
  nbeads = src.nbeads;
  nstates = src.nstates;
  nnucl = src.nnucl;

  ovlp_dia = src.ovlp_dia.clone();
  ham_dia = src.ham_dia.clone();
  nac_dia = src.nac_dia.clone();
  hvib_dia = src.hvib_dia.clone();
  dc1_dia = src.dc1_dia.clone();
  d1ham_dia = src.d1ham_dia.clone();
  d2ham_dia = src.d2ham_dia.clone();

  ovlp_adi = src.ovlp_adi.clone();
  ham_adi = src.ham_adi.clone();
  nac_adi = src.nac_adi.clone();
  hvib_adi = src.hvib_adi.clone();
  dc1_adi = src.dc1_adi.clone();
  d1ham_adi = src.d1ham_adi.clone();
  d2ham_adi = src.d2ham_adi.clone();

  basis_transform = src.basis_transform.clone();
  time_overlap_dia = src.time_overlap_dia.clone();
  time_overlap_adi = src.time_overlap_adi.clone();

}


sHamiltonian::sHamiltonian(int nbeads_, int nelec_, int nnucl_){

  nbeads = nbeads_;
  nstates = nelec_; 
  nnucl = nnucl_;
 
  ovlp_dia = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  ham_dia = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  nac_dia = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  hvib_dia = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  dc1_dia = torch::zeros({ nbeads, nnucl, nstates, nstates }, torch::kComplexDouble);
  d1ham_dia = torch::zeros({ nbeads, nnucl, nstates, nstates }, torch::kComplexDouble);
  d2ham_dia = torch::zeros({ nbeads, nnucl, nnucl, nstates, nstates }, torch::kComplexDouble);

  ovlp_adi = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  ham_adi = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  nac_adi = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  hvib_adi = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  dc1_adi = torch::zeros({ nbeads, nnucl, nstates, nstates }, torch::kComplexDouble);
  d1ham_adi = torch::zeros({ nbeads, nnucl, nstates, nstates }, torch::kComplexDouble);
  d2ham_adi = torch::zeros({ nbeads, nnucl, nnucl, nstates, nstates }, torch::kComplexDouble);

  basis_transform = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  time_overlap_dia = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);
  time_overlap_adi = torch::zeros({ nbeads, nstates, nstates }, torch::kComplexDouble);

}



sHamiltonian::~sHamiltonian(){ 
  ;;
}


void sHamiltonian::bind(std::string name, torch::Tensor x){

  if(name=="ovlp_dia"){  ovlp_dia = x; }
  else if(name=="ham_dia"){  ham_dia = x; }
  else if(name=="nac_dia"){  nac_dia = x; }
  else if(name=="hvib_dia"){  hvib_dia = x; }
  else if(name=="dc1_dia"){  dc1_dia = x; }
  else if(name=="d1ham_dia"){  d1ham_dia = x; }
  else if(name=="d2ham_dia"){  d2ham_dia = x; }

  else if(name=="ovlp_adi"){  ovlp_adi = x; }
  else if(name=="ham_adi"){  ham_adi = x; }
  else if(name=="nac_adi"){  nac_adi = x; }
  else if(name=="hvib_adi"){  hvib_adi = x; }
  else if(name=="dc1_adi"){  dc1_adi = x; }
  else if(name=="d1ham_adi"){  d1ham_adi = x; }
  else if(name=="d2ham_adi"){  d2ham_adi = x; }

  else if(name=="basis_tansform"){  basis_transform = x; }
  else if(name=="time_overlap_dia"){  time_overlap_dia = x; }
  else if(name=="time_overlap_adi"){  time_overlap_adi = x; }

}



void sHamiltonian::compute(std::string property, bp::object py_funct, torch::Tensor q, bp::object params){

  if(property=="ovlp_dia"){    ovlp_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="ham_dia"){  ham_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="nac_dia"){    nac_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="hvib_dia"){    hvib_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="dc1_dia"){    dc1_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="d1ham_dia"){    d1ham_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="d2ham_dia"){    d2ham_dia = extract<torch::Tensor>(py_funct(q, params));    }

  else if(property=="ovlp_adi"){    ovlp_adi = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="ham_adi"){  ham_adi = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="nac_adi"){    nac_adi = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="hvib_adi"){    hvib_adi = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="dc1_adi"){    dc1_adi = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="d1ham_adi"){    d1ham_adi = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="d2ham_adi"){    d2ham_adi = extract<torch::Tensor>(py_funct(q, params));    }

  else if(property=="basis_transform"){    basis_transform = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="time_overlap_dia"){   time_overlap_dia = extract<torch::Tensor>(py_funct(q, params));    }
  else if(property=="time_overlap_adi"){   time_overlap_adi = extract<torch::Tensor>(py_funct(q, params));    }

}


void sHamiltonian::compute(bp::object py_funct, torch::Tensor q, bp::object params){

  bp::object obj = py_funct(q, params);


  if(  (int)hasattr(obj,"ovlp_dia") ){    ovlp_dia =  extract<torch::Tensor>( obj.attr("ovlp_dia"));    }
  if(  (int)hasattr(obj,"ham_dia") ){    ham_dia =  extract<torch::Tensor>( obj.attr("ham_dia"));    }
  if(  (int)hasattr(obj,"nac_dia") ){    nac_dia =  extract<torch::Tensor>( obj.attr("nac_dia"));    }
  if(  (int)hasattr(obj,"hvib_dia") ){    hvib_dia =  extract<torch::Tensor>( obj.attr("hvib_dia"));    }
  if(  (int)hasattr(obj,"dc1_dia") ){    dc1_dia =  extract<torch::Tensor>( obj.attr("dc1_dia"));    }
  if(  (int)hasattr(obj,"d1ham_dia") ){    d1ham_dia =  extract<torch::Tensor>( obj.attr("d1ham_dia"));    }
  if(  (int)hasattr(obj,"d2ham_dia") ){    d2ham_dia =  extract<torch::Tensor>( obj.attr("d2ham_dia"));    }

  if(  (int)hasattr(obj,"ovlp_adi") ){    ovlp_adi =  extract<torch::Tensor>( obj.attr("ovlp_adi"));    }
  if(  (int)hasattr(obj,"ham_adi") ){    ham_adi =  extract<torch::Tensor>( obj.attr("ham_adi"));    }
  if(  (int)hasattr(obj,"nac_adi") ){    nac_adi =  extract<torch::Tensor>( obj.attr("nac_adi"));    }
  if(  (int)hasattr(obj,"hvib_adi") ){    hvib_adi =  extract<torch::Tensor>( obj.attr("hvib_adi"));    }
  if(  (int)hasattr(obj,"dc1_adi") ){    dc1_adi =  extract<torch::Tensor>( obj.attr("dc1_adi"));    }
  if(  (int)hasattr(obj,"d1ham_adi") ){    d1ham_adi =  extract<torch::Tensor>( obj.attr("d1ham_adi"));    }
  if(  (int)hasattr(obj,"d2ham_adi") ){    d2ham_adi =  extract<torch::Tensor>( obj.attr("d2ham_adi"));    }

  if(  (int)hasattr(obj,"basis_transform") ){    basis_transform =  extract<torch::Tensor>( obj.attr("basis_transform"));  }
  if(  (int)hasattr(obj,"time_overlap_dia") ){    time_overlap_dia =  extract<torch::Tensor>( obj.attr("time_overlap_dia"));    }
  if(  (int)hasattr(obj,"time_overlap_adi") ){    time_overlap_adi =  extract<torch::Tensor>( obj.attr("time_overlap_adi"));    }


}


void sHamiltonian::dia2adi(){

    const int64_t batch_size = ham_dia.size(0);
    const int64_t n = ham_dia.size(1);


// Compute eigenvalues and eigenvectors
    //auto result = torch::linalg_eig(ham_dia);
    //auto eigenvalues = std::get<0>(result);
    //auto eigenvectors = std::get<1>(result);

    // Extract real part of eigenvalues
    //auto eigenvalues_real = torch::real(eigenvalues);

    // Sort eigenvalues and get sorted indices
    //torch::Tensor sorted_eigenvalues;
    //torch::Tensor sorted_indices;
    //std::tie(sorted_eigenvalues, sorted_indices) = torch::sort(eigenvalues_real, 1);

    // Reorder eigenvalues and eigenvectors
    //sorted_eigenvalues = torch::gather(eigenvalues, 1, sorted_indices);
    //auto sorted_eigenvectors = torch::gather(eigenvectors, 2, sorted_indices.unsqueeze(-1).expand({batch_size, n, n}));


    
    // 1. Compute batched eigendecomposition
    auto [eigenvalues, eigenvectors] = torch::linalg_eig(ham_dia); 
    // eigenvalues shape: [batch_size, N] (complex)
    // eigenvectors shape: [batch_size, N, N] (complex)

    // 2. Compute sorting keys (magnitudes)
    //auto magnitudes = torch::abs(eigenvalues); // [batch_size, N]
    auto eig_re = torch::real(eigenvalues);
    auto eig_im = torch::imag(eigenvalues);



    // 3. Get sorted indices and eigenvectors
    auto [sorted_eigenvalues_re, sorted_indices] = torch::sort(eig_re, 1, false);

    // Sort the imag part of the eigenvalues:
    auto batch_idx = torch::arange(batch_size).unsqueeze(1);
    auto sorted_eigenvalues_im = torch::zeros_like(sorted_eigenvalues_re);
    sorted_eigenvalues_im.index_put_( {batch_idx, sorted_indices},  eig_im );

    auto sorted_eigenvalues = torch::complex(sorted_eigenvalues_re, sorted_eigenvalues_im);
    
    // 4. Gather the sorted eigenvectors
    // indices: (B, N) -> expand to (B, 1, N) for batch-wise gathering
    auto indices_expanded = sorted_indices.unsqueeze(1).expand({-1, eigenvectors.size(1), -1});

    // Gather eigenvectors along the last dimension
    auto sorted_eigenvectors = eigenvectors.gather(2, indices_expanded);


/*
  cout<<"============== Real components ===============\n";
  cout<<"==== Eigenvalues ======\n";
  cout<<"eigenvalues = "<<torch::real(eigenvalues)<<endl;
  cout<<"sorted eigenvalues = "<<sorted_eigenvalues_re<<endl;

  cout<<"==== Eigenvectors ======\n";
  cout<<"eigenvectors = "<<torch::real(eigenvectors)<<endl;
  cout<<"sorted eigenvalues = "<<torch::real(sorted_eigenvectors)<<endl;

  cout<<"============== Imag components ===============\n";
  cout<<"==== Eigenvalues ======\n";
  cout<<"eigenvalues = "<<torch::imag(eigenvalues)<<endl;
  cout<<"sorted eigenvalues = "<<sorted_eigenvalues_im<<endl;

  cout<<"==== Eigenvectors ======\n";
  cout<<"eigenvectors = "<<torch::imag(eigenvectors)<<endl;
  cout<<"sorted eigenvectors = "<<torch::imag(sorted_eigenvectors)<<endl;
*/

  auto ham_adi_real = torch::diag_embed( sorted_eigenvalues_re );
  auto ham_adi_imag = torch::diag_embed( sorted_eigenvalues_im );  
  

  ham_adi = torch::complex(ham_adi_real, ham_adi_imag);
  basis_transform = sorted_eigenvectors;


 
}


void sHamiltonian::compute_nacs_and_grads(){

  // E.g. see the derivations here: https://github.com/alexvakimov/Derivatory/blob/master/theory_NAC.pdf^M
  // also: http://www.theochem.ruhr-uni-bochum.de/~nikos.doltsinis/nic_10_doltsinis.pdf^M
  auto tmp = torch::matmul( basis_transform.conj().transpose(-2, -1), torch::matmul(d1ham_dia, basis_transform) );

//        *dtilda = (*basis_transform) * (*dc1_dia[n]).H() * (*basis_transform).H() * (*ham_adi);^M
  auto dtilda = torch::matmul( basis_transform.conj().transpose(-2, -1), torch::matmul(dc1_dia, torch::matmul( basis_transform, ham_adi)));
  dtilda = dtilda + dtilda.conj().transpose(-2, -1);
  tmp = tmp - dtilda;

  // Adiabatic "forces"
  d1ham_adi.zero_();
  dc1_adi.zero_();

  //auto ham_acc = ham_adi.accessor< torch::kComplexDouble, 3>();  // 3D complex tensor

  for(int b = 0; b < nbeads; b++){
    for(int n = 0; n < nnucl; n++){
      for(int i = 0; i<nstates; i++){  
 
        d1ham_adi[b][n][i][i] = tmp[b][n][i][i];

        for(int j = i+1; j < nstates; j++){
          double dE = at::real(ham_adi[b][j][j] - ham_adi[b][i][i]).item<double>();
          if(fabs(dE)<1e-25){ dE = 1e-25; }

          dc1_adi[b][n][i][j] = tmp[b][n][i][j]/dE;
          dc1_adi[b][n][j][i] = -dc1_adi[b][n][i][j];
        }// for j
      }// for i
    }// for n
  }// b

}// compute_nacs_and_grads

torch::Tensor sHamiltonian::forces_adi(){
// [nbeads, nstates, nnucl]
  torch::Tensor f = torch::zeros({ nbeads, nstates, nnucl }, torch::kDouble); 

  for(int b = 0; b < nbeads; b++){
    for(int n = 0; n < nnucl; n++){
      for(int i = 0; i<nstates; i++){
        f[b][i][n] = -at::real(d1ham_adi[b][n][i][i]).item<double>();
      }// for i
    }// for n
  }// for b

  return f;
}

torch::Tensor sHamiltonian::forces_dia(){
// [nbeads, nstates, nnucl]
  torch::Tensor f = torch::zeros({ nbeads, nstates, nnucl }, torch::kDouble);

  for(int b = 0; b < nbeads; b++){
    for(int n = 0; n < nnucl; n++){
      for(int i = 0; i<nstates; i++){
        f[b][i][n] = -at::real(d1ham_dia[b][n][i][i]).item<double>();
      }// for i
    }// for n
  }// for b

  return f;
}




}// namespace libshamiltonian
}// liblibra

