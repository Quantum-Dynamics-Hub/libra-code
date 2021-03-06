#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *




class patt():
    pass
p1 = patt()
p1.Input = [0.0, 0.0]
p1.Output = [0.0]
p1.Derivs = [ 0.0,0.0 ]
print p1

p2 = patt()
p2.Input = [0.5,0.0]
p2.Output = [0.0]
p2.Derivs = [ 0.0,0.0 ]
print p2

p3 = patt()
p3.Input = [0.0,0.5]
p3.Output = [0.0]
p3.Derivs = [ 0.0,0.0 ]
print p3

p4 = patt()
p4.Input = [0.5,0.5]
p4.Output = [0.5]
p4.Derivs = [ 0.0,0.0 ]
print p4

training_set = [ p1,p4] # p2,p3,p4 ]


num_of_patterns = 4
epoch_size = 4
rperm = intList()
for i in xrange(10):
    randperm(epoch_size,num_of_patterns,rperm)
    print Cpp2Py(rperm)

#sys.exit(0)

print "\nTest 1: Running setup from other test..."
ANN = NeuralNetwork()
ANN.CreateANN([2, 2, 2, 1])


#print ANN.B
#print ANN.Inputs
#print ANN.Outputs

#sys.exit()

ANN.SetTrainingData(training_set,1)  # 0 - means we don't really use info about derivatives, 1 - we use it, but we have to use 1
#ANN.NormalizeAndScaleTrainingData(1,1,[-0.5, -0.5, -0.5], [0.5, 0.5, 0.5])
#ANN.NormalizeAndTransformTrainingData(1,1)
ANN.ScaleTrainingData(1,1)

#sys.exit()

print "\nTest2: Now set configuration"
print "class ann_config_class():"
print "    pass"
print "ann_config = ann_config_class()"
print "ann_config.learning_method = \"RProp\"  # other options: BackProp, QuickProp, ConjGradProp"
print "ANN.set(ann_config)"


class ann_config_class():
    pass
ann_config = ann_config_class()
ann_config.learning_method = "BackProp"  # all options: RProp, BackProp, QuickProp, ConjGradProp
ann_config.learning_rate = 0.01
ann_config.epoch_size = 2
ann_config.momentum_term = 0.0
ann_config.grad_weight = 0.0
ann_config.norm_exp = 0.5
ann_config.iterations_in_cycle = 1000
ANN.set(ann_config)

#sys.exit()

print "\nTest3: Start training process"
print" ANN.ANNTrain()"
ANN.ANNTrain()

print "\nTest4: Learning history"
print "The output is organized like this: there n_inp lists - one for each expected output, then we have"
print "n_inp x n_out lists - one for each derivatiave of output w.r.t. inputs"
print "Each list contains 5 elements: the instantaneous error for given channel, the average error [percentage of training output], max error [percentage of training output], minimal error [percentage], percentage of the patterns that are approximated with less than 50% error"
ANN.LearningHistory("hist.txt","original")

print "\nTest5: Continue training in loop"
for i in range(0,50):
    ANN.ANNTrain()
    ANN.LearningHistory("hist1.txt","original")
    ANN.ExportANN("ann"+str(i))    


print "\nTest6: Predictive capabilities of the trained ANN"
for p in training_set:
    res = []
    ANN.Propagate(p.Input, res)
    print p.Input, " --> ", "Actual= ", res, " Expected= ", p.Output



