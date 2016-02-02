/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

//#include <string>
#include <boost/python.hpp>
#include "libann.h"

using namespace boost::python;

/// libmmath namespace
namespace libmmath{

/// libann namespace
namespace libann{


void export_NeuralNetwork_objects(){



int (NeuralNetwork::*ScaleTrainingData1)(int,int)                                         = &NeuralNetwork::ScaleTrainingData;
int (NeuralNetwork::*ScaleTrainingData2)(int,int,boost::python::list,boost::python::list) = &NeuralNetwork::ScaleTrainingData;

int (NeuralNetwork::*Propagate1)(MATRIX,MATRIX&)                                  = &NeuralNetwork::Propagate;
int (NeuralNetwork::*Propagate2)(boost::python::list,boost::python::list&)        = &NeuralNetwork::Propagate;
int (NeuralNetwork::*Propagate3)(MATRIX,MATRIX&,MATRIX&)                          = &NeuralNetwork::Propagate;
int (NeuralNetwork::*Propagate4)(boost::python::list,boost::python::list&,boost::python::list&)        = &NeuralNetwork::Propagate;


//---------------------------------------------------------------------



//-------- Write documentation of all classes, functions and class members here ------------------------
std::string tmp;

std::string set_docstring;
std::string set_training_data_docstring;
std::string NormalizeTrainingData_docstring;
std::string ScaleTrainingData1_docstring;
std::string ScaleTrainingData2_docstring;
std::string ImportANN_docstring;
std::string ExportANN_docstring;


set_docstring="";       tmp = "Creates an ANN from the corresponding architechture object\n";
set_docstring += tmp;   tmp = "Function prototype set(object ann)\n";
set_docstring += tmp;   tmp = " ann - is a python class defining the ANN architechture. It should \n";
set_docstring += tmp;   tmp = " contain 3 data fields: Title, NPE, Bias (python lists)\n";
set_docstring += tmp;   tmp = " Title - definition of each corresponding entry of NPE and Bias lists\n";
set_docstring += tmp;   tmp = " NPE - the number of neurons in each layer(the number of layers is a list size)\n";
set_docstring += tmp;   tmp = "       it includes input neurons and output neurons\n";
set_docstring += tmp;   tmp = " Bias - originally was included in order to initialize biases, but may be omitted\n";
set_docstring += tmp;   tmp = " Example:\n";
set_docstring += tmp;   tmp = " class aux():\n";
set_docstring += tmp;   tmp = "     pass\n";
set_docstring += tmp;   tmp = " ann = aux() \n";
set_docstring += tmp;   tmp = " ANN = NeuralNetwork()\n";
set_docstring += tmp;   tmp = " ann.Title = ['Input','First hidden layer','Second hidden layer','Output']\n";
set_docstring += tmp;   tmp = " ann.NPE   = [3,5,5,4]\n";
set_docstring += tmp;   tmp = " ANN.set(ann)\n";
set_docstring += tmp;   tmp = " This will create an ANN with 2 hidden layers with 5 neurons in each of the hidden layers\n";
set_docstring += tmp;   tmp = " and 3 input and 4 output neurons\n";
set_docstring += tmp; 

set_training_data_docstring="";       tmp = "Function which presents the training data to the ANN\n";
set_training_data_docstring += tmp;   tmp = "Function prototype set(object obj)\n";
set_training_data_docstring += tmp;   tmp = "  Need to set the training_set object as follows:\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = "  obj - is a list of patterns, where each pattern is the\n";
set_training_data_docstring += tmp;   tmp = "        list of two other lists - one of them - is input\n";
set_training_data_docstring += tmp;   tmp = "        another list - is output\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = "  obj = [ pattern1,\n";
set_training_data_docstring += tmp;   tmp = "          pattern2,\n";
set_training_data_docstring += tmp;   tmp = "            ...\n";
set_training_data_docstring += tmp;   tmp = "          patternN\n";
set_training_data_docstring += tmp;   tmp = "        ]\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = "  pattern.Input  = [ vector of size n]\n";
set_training_data_docstring += tmp;   tmp = "  pattern.Output = [ vector of size m]\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = "  This is for the ANN architecture:\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = " Input                              Output\n";
set_training_data_docstring += tmp;   tmp = " layer                              layer\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = "    1                                  1\n";
set_training_data_docstring += tmp;   tmp = "    2                                  2\n";
set_training_data_docstring += tmp;   tmp = "   ...   --->   Hidden layers --->    ...\n";
set_training_data_docstring += tmp;   tmp = "\n";
set_training_data_docstring += tmp;   tmp = "    n                                  m\n";
set_training_data_docstring += tmp;   tmp = "Example:\n";
set_training_data_docstring += tmp;   tmp = "    pattern = aux()\n";
set_training_data_docstring += tmp;   tmp = "    pattern.Input = [1,1]\n";
set_training_data_docstring += tmp;   tmp = "    pattern.Output = [1]\n";
set_training_data_docstring += tmp;   tmp = "    training_set.append(pattern)\n";
set_training_data_docstring += tmp;   tmp = "    ....\n";
set_training_data_docstring += tmp;   tmp = "    ANN.set_training_data(training_set)\n";
set_training_data_docstring += tmp;   tmp = " This will add 1(or more) patterns, where each pattern has 2 input data\n";
set_training_data_docstring += tmp;   tmp = " and 1 output data\n";
set_training_data_docstring += tmp;

NormalizeTrainingData_docstring ="";     tmp = "Performs linear scaling of both input and output vectors according to formula:\n";
NormalizeTrainingData_docstring += tmp;  tmp = " x[i] <- (x[i]-<x>)/sigma \n";
NormalizeTrainingData_docstring += tmp;  tmp = " where <x> - is an average value of variable x: <x> = (1/N)summ_over_i(x[i])\n";
NormalizeTrainingData_docstring += tmp;  tmp = " and sigma - is a variance: sigma^2 = (1/(N-1))summ_over_i((x[i]-<x>)^2)\n";
NormalizeTrainingData_docstring += tmp;  tmp = " Arguments: none\n";
NormalizeTrainingData_docstring += tmp;  tmp = " Note that this kind of regularization does not necessarily leads to target values in range of [-1,1]\n";
NormalizeTrainingData_docstring += tmp; 

ScaleTrainingData1_docstring ="";     tmp = "Performs linear scaling of both input and output vectors depending of X_scale and Y_scale flags\n";
ScaleTrainingData1_docstring += tmp;  tmp = "This will scale corresponding vectors onto the [-0.95,0.95] interval\n";
ScaleTrainingData1_docstring += tmp;  tmp = "Arguments: none\n";
ScaleTrainingData1_docstring += tmp;

ScaleTrainingData2_docstring ="";    tmp = "Performs linear scaling of both input and output vectors\n";
ScaleTrainingData2_docstring += tmp; tmp = "This will scale corresponding vectors onto the intervals defined by upper and lower limits \n";
ScaleTrainingData2_docstring += tmp; tmp = "defined by minlim and maxlim lists\n";
ScaleTrainingData2_docstring += tmp; tmp = "Arguments: minlim - minimal boundary of the interval to which we want to map original data\n";
ScaleTrainingData2_docstring += tmp; tmp = "           maxlim - maximal boundary of the interval to which we want to map original data\n";
ScaleTrainingData2_docstring += tmp; tmp = "If minimal and maximal boundaries are zero - the will be no scaling of corresponding component of the pattern\n";
ScaleTrainingData2_docstring += tmp; tmp = "Example:\n";
ScaleTrainingData2_docstring += tmp; tmp = "ANN.ScaleTrainingData([-5.0,-3.0,-1.0],[5.0,6.0,1.0])\n";
ScaleTrainingData2_docstring += tmp; tmp = "is valid for patterns of 2 inputs and 1 output (or 1 input or 2 outputs) and will\n";
ScaleTrainingData2_docstring += tmp; tmp = "scale all 1-st components of the pattern vectors onto interval [-5.0,5.0]\n";
ScaleTrainingData2_docstring += tmp; tmp = "scale all 2-nd components of the pattern vectors onto interval [-3.0,6.0]\n";
ScaleTrainingData2_docstring += tmp; tmp = "scale all 3-rd components of the pattern vectors onto interval [-1.0,1.0]\n";
ScaleTrainingData2_docstring += tmp; tmp = "Note, however, the intervals to which output components of the patterns will be mapped should be\n";
ScaleTrainingData2_docstring += tmp; tmp = "chosen according to the range of the ANN output, that is [-1.0,1.0]\n";
ScaleTrainingData2_docstring += tmp;

ImportANN_docstring = "";            tmp = "Creates or updates the neural network using a file of specific format and given as argument\n";
ImportANN_docstring += tmp;          tmp = "If the object is new - this will create a new fully-functional neural network\n";
ImportANN_docstring += tmp;          tmp = "If the object existed - old information will be erased and new one will be added\n";
ImportANN_docstring += tmp;          tmp = "Arguments: filename - the name of the file from which to read the information\n";
ImportANN_docstring += tmp;          tmp = "           the format of the 'filename' is consistent with the format of ExportANN function\n";
ImportANN_docstring += tmp;        

ExportANN_docstring = "";         tmp = "Writes the current state of the neural network to the file of specific format\n";
ExportANN_docstring += tmp;       tmp = "The information writen includes(in this order): number of layers, network architecture\n";
ExportANN_docstring += tmp;       tmp = "weights,biases,scaling parameters and the scaling method used\n";
ExportANN_docstring += tmp;       tmp = "Arguments: filename - the name of the file to which write the information\n";
ExportANN_docstring += tmp;       tmp = "           the format of the 'filename' made recognizable by ImportANN function\n";
ExportANN_docstring += tmp;       tmp = "           so the resulting file may be further used for fast and easy construction of the ANN\n";
ExportANN_docstring += tmp;       tmp = "           or to continue training from the best state obtained (not from the random guess)\n";
ExportANN_docstring += tmp;




//---------------------------------------------------------------------

    class_<NeuralNetwork>("NeuralNetwork",init<>())
        .def("CreateANN",&NeuralNetwork::CreateANN)
        .def("ShowANN",&NeuralNetwork::ShowANN)


        .def("ExportANN",&NeuralNetwork::ExportANN,ExportANN_docstring.c_str())
        .def("ImportANN",&NeuralNetwork::ImportANN,ImportANN_docstring.c_str())
        .def("set",&NeuralNetwork::set,set_docstring.c_str())
     
        .def("SetTrainingData",&NeuralNetwork::SetTrainingData,set_training_data_docstring.c_str())      
        .def("NormalizeTrainingData",&NeuralNetwork::NormalizeTrainingData,NormalizeTrainingData_docstring.c_str())
        .def("ScaleTrainingData",ScaleTrainingData1,ScaleTrainingData1_docstring.c_str())
        .def("ScaleTrainingData",ScaleTrainingData2,ScaleTrainingData2_docstring.c_str())
        .def("NormalizeAndScaleTrainingData",&NeuralNetwork::NormalizeAndScaleTrainingData)
        .def("NormalizeAndTransformTrainingData", &NeuralNetwork::NormalizeAndTransformTrainingData)
        .def("CropTrainingData",&NeuralNetwork::CropTrainingData)

        .def("Propagate",Propagate1)
        .def("Propagate",Propagate2)
        .def("Propagate",Propagate3)
        .def("Propagate",Propagate4)
        .def("ANNTrain",&NeuralNetwork::ANNTrain)
        .def("LearningHistory",&NeuralNetwork::LearningHistory)


        .def_readwrite("Inputs",&NeuralNetwork::Inputs)
        .def_readwrite("Outputs",&NeuralNetwork::Outputs)
        .def_readwrite("B",&NeuralNetwork::B)
         
    ;

}




#ifdef CYGWIN
BOOST_PYTHON_MODULE(cygann){
#else
BOOST_PYTHON_MODULE(libann){
#endif

  //to_python_converter<std::vector<libmmath::DATA>, VecToList<libmmath::DATA> >();

//  export_Mathematics_objects();  // also register mmath python functions!
  export_NeuralNetwork_objects();

}


}// namespace libann
}// namespace libmmath




