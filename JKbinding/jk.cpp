#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>
#include<omp.h>
#include<sstream>
#include<stdio.h>
#include<vector>
#include<math.h>


namespace py = pybind11;

int jk_numpy(py::array_t<double> I, py::array_t<double> D,
             py::array_t<double> J, py::array_t<double> K)
{
  py::buffer_info I_info = I.request();
  py::buffer_info D_info = D.request();
  py::buffer_info J_info = J.request();
  py::buffer_info K_info = K.request();

  size_t d = I_info.shape[0];
  size_t d2 = d * d;
  size_t d3 = d * d * d;

  const double * I_data = static_cast<double *>(I_info.ptr);
  const double * D_data = static_cast<double *>(D_info.ptr);
  double * J_data = static_cast<double *>(J_info.ptr);
  double * K_data = static_cast<double *>(K_info.ptr);

#pragma omp parallel for schedule(dynamic) num_threads(2)  //must specify the number of threads here
for(size_t p = 0; p < d; p++)
  {
    for(size_t q = 0; q <= p; q++)
    {
      double valj = 0.0;
      double valk = 0.0;
      int ind = p * d3 + q * d2, ind2 = p * d3 + q * d;
      for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j <= i; j++) {
          valj += 2.0 * I_data[ind + i* d +j]  * D_data[i * d + j];
          valk += (I_data[ind2 + j * d2 + i] + I_data[ind2 + i * d2 + j]) * D_data[i * d + j];
        }
      }
      J_data[p * d + q] = valj;
      J_data[q * d + p] = valj;
      K_data[p * d + q] = valk;
      K_data[q * d + p] = valk;
    }
  }

  return 0;
}

PYBIND11_PLUGIN(qm7_jk)
{
  py::module m("qm7_jk", "QM 7 JK");

  m.def("jk_numpy", &jk_numpy, "JK");

  return m.ptr();
}
