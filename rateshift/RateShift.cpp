//
// File: RateShift.cpp
// Created by: Julien Dutheil
// Created on: Thrusday, February 27th 2020 22:19
//

/*
    This file is part of the BppRateShift package.

    DiffRate is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PhySamp is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

#include <Bpp/App/BppApplication.h>

// From bpp-seq:
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>
#include <Bpp/Phyl/Legacy/Likelihood/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bpprateshift parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to http://BioPP/rateshift for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*               Bio++ Rate Shift, version 1.0.0.                 *" << endl;
  cout << "* Author: J. Y. Dutheil                     Last Modif. 04/09/23 *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Original work: Tal Pupko & Nicolas Galtier                     *" << endl;
  cout << "*                     Proc Biol Sci. 2002 Jul 7;269(1498):1313-6.*" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication rateshift(args, argv, "DiffRate");
  rateshift.startTimer();

  //Get alphabet and genetic code (if needed):
  shared_ptr<const Alphabet> alphabet = SequenceApplicationTools::getAlphabet(rateshift.getParams());
  shared_ptr<const GeneticCode> gCode;
  auto codonAlphabet = dynamic_pointer_cast<const CodonAlphabet>(alphabet);
  if (codonAlphabet) {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", rateshift.getParams(), "Standard", "", true, true);
    ApplicationTools::displayResult("Genetic Code", codeDesc);
      
    gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
  }

  //Get alignment:
  shared_ptr<SiteContainerInterface> sites = SequenceApplicationTools::getSiteContainer(alphabet, rateshift.getParams());

  //Get tree:
  auto tmpTree = PhylogeneticsApplicationTools::getTree(rateshift.getParams());
  auto tree = make_unique<TreeTemplate<Node>>(*tmpTree);

  //Eventually, print tree with id to a file in order to select foreground branches:
  string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", rateshift.getParams(), false, false, "", true, "none", 1);
  if (treeWIdPath != "none")
  {
    vector<Node*> nodes = tree->getNodes();
    for (auto node : nodes)
    {
      if (node->isLeaf())
        node->setName(TextTools::toString(node->getId()) + "_" + node->getName());
      else
        node->setBranchProperty("NodeId", BppString(TextTools::toString(node->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.writeTree(*tree, treeWIdPath, true);
    cout << "BppRateShift's done." << endl;
    exit(0);
  }



  //Get substitution model:
  shared_ptr<const AlignmentDataInterface> data = sites;
  map<string, string> unparsedParams;
  shared_ptr<TransitionModelInterface> model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode, data, rateshift.getParams(), unparsedParams);
  shared_ptr<DiscreteDistribution> rDist = PhylogeneticsApplicationTools::getRateDistribution(rateshift.getParams());

  //Optional: optimize model parameters first
  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
  auto htl = shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface>(new RHomogeneousTreeLikelihood(*tree, *data, model, rDist, true, true, false));
  htl->initialize();

  htl = dynamic_pointer_cast<DiscreteRatesAcrossSitesTreeLikelihoodInterface>(
    PhylogeneticsApplicationToolsOld::optimizeParameters(htl, htl->getParameters(), rateshift.getParams()));

  tree.reset(new TreeTemplate<Node>(htl->tree()));
  PhylogeneticsApplicationTools::writeTree(*tree, rateshift.getParams());

  // Write parameters to screen:
  ApplicationTools::displayResult("Log likelihood", TextTools::toString(htl->getValue(), 15));
  ParameterList parameters = htl->getSubstitutionModelParameters();
  for (size_t i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  parameters = htl->getRateDistributionParameters();
  for (size_t i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }

  // Checking convergence:
  PhylogeneticsApplicationTools::checkEstimatedParameters(htl->getParameters());

  //Add rate parameter:
  model->addRateParameter();
  string rateParam = model->getNamespace() + "rate";

  //Use a constant rate:
  rDist.reset(new ConstantRateDistribution());

  //Create a one-rate tree likelihood:
  shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> oneRateTl(
      new RHomogeneousTreeLikelihood(*tree, model, rDist, true, false, false));

  //Create a two-rate tree likelihood:
  //First, create a model set with two models, according to the foreground branches:
  vector<int> foregroundIds = ApplicationTools::getVectorParameter<int>("foreground_branches", rateshift.getParams(), ',', '-', "");
  vector<int> backgroundIds;
  vector<int> allIds = tree->getBranchesId();
  for (auto id: allIds) {
    if (! VectorTools::contains(foregroundIds, id))
      backgroundIds.push_back(id);
  }
  ApplicationTools::displayResult("Number of foreground branches", foregroundIds.size());
  ApplicationTools::displayResult("Number of background branches", backgroundIds.size());

  shared_ptr<SubstitutionModelSet> modelSet(new SubstitutionModelSet(alphabet));
  modelSet->addModel(shared_ptr<TransitionModelInterface>(model->clone()), foregroundIds);
  modelSet->addModel(shared_ptr<TransitionModelInterface>(model->clone()), backgroundIds);
  string rateParam1 = model->getNamespace() + "rate_1";
  string rateParam2 = model->getNamespace() + "rate_2";
  vector<string> rateParams{rateParam1, rateParam2};
  shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> twoRateTl(
      new RNonHomogeneousTreeLikelihood(*tree, modelSet, rDist, false, false, false));
  

  //Loop over all sites and perform a likelihood ratio test. Print to SGED file.
  string output = ApplicationTools::getAFilePath("output.file", rateshift.getParams(), true, false);
  ApplicationTools::displayResult("Writing results to", output);
  ofstream out(output, ios::out);
  out << "Group\tr\tr.fg\tr.bg\tAIC1\tAIC2\tdiffLnL\tP.value" << endl;
  ApplicationTools::displayTask("Testing all sites", true);
  for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
    ApplicationTools::displayGauge(i, sites->getNumberOfSites() - 1, '=');
    SiteSelection s;
    s.push_back(i);
    auto site = SiteContainerTools::getSelectedSites(*sites, s);

    oneRateTl->setData(*site);
    oneRateTl->initialize();
    OptimizationToolsOld::optimizeNumericalParameters2(oneRateTl,
        oneRateTl->getParameters().createSubList(rateParam),
        0, 0.000001, 10000, nullptr, nullptr, false, false, 0,
        OptimizationTools::OPTIMIZATION_NEWTON);
    
    twoRateTl->setData(*site);
    twoRateTl->initialize();
    OptimizationToolsOld::optimizeNumericalParameters2(twoRateTl,
        twoRateTl->getParameters().createSubList(rateParams),
        0, 0.000001, 10000, nullptr, nullptr, false, false, 0,
        OptimizationTools::OPTIMIZATION_NEWTON);

    double diffLnL = twoRateTl->getLogLikelihood() - oneRateTl->getLogLikelihood();
    double pvalue = 1. - RandomTools::pChisq(2*diffLnL, 1.);
    double aic1 = 2. - 2.*oneRateTl->getLogLikelihood();
    double aic2 = 4. - 2.*twoRateTl->getLogLikelihood();


    out << "[" << sites->site(i).getCoordinate() << "]" << "\t";
    out << model->getParameters().getParameterValue(rateParam);
    out << "\t";
    out << modelSet->getParameters().getParameterValue(rateParam1);
    out << "\t";
    out << modelSet->getParameters().getParameterValue(rateParam2);
    out << "\t";
    out << aic1;
    out << "\t";
    out << aic2;
    out << "\t";
    out << diffLnL;
    out << "\t";
    out << pvalue;
    out << endl;
  }
  ApplicationTools::displayTaskDone();

  //Exiting...
  rateshift.done();
  }
  catch (exception& e)
  {
    cout << endl;
    cout << "_____________________________________________________" << endl;
    cout << "ERROR!!!" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

